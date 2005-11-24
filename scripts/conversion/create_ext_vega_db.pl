#!/usr/local/bin/perl

=head1 NAME

create_ext_vega_db.pl - filter internal annotation from a Vega db

=head1 SYNOPSIS

create_ext_vega_db.pl [options]

General options:
    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)

    --dbname, db_name=NAME              use database NAME
    --host, --dbhost, --db_host=HOST    use database host HOST
    --port, --dbport, --db_port=PORT    use database port PORT
    --user, --dbuser, --db_user=USER    use database username USER
    --pass, --dbpass, --db_pass=PASS    use database passwort PASS
    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)
    -v, --verbose                       verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)

Specific options:

    --infile=FILE                       read list of stable IDs from FILE
    
=head1 DESCRIPTION


=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options('infile=s');
$support->allowed_params($support->get_common_params, qw(infile));

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

$support->check_required_params('infile');

# connect to database and get adaptors
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;

# read list of stable IDs to keep in external db - everything else will be
# deleted
$support->log_stamped("Reading stable IDs to keep from file...\n");
my $in = $support->filehandle('<', $support->param('infile'));
my @gene_stable_ids;
while (<$in>) {
    chomp;
    push @gene_stable_ids, $_;
}
close($in);
my $gsi_list = "'" . join("', '", @gene_stable_ids) . "'";
$support->log_stamped("Done reading ".scalar(@gene_stable_ids)." entries.\n\n");

if ($support->param('dry_run')) {
    $support->log("Nothing else to be done for dry run. Aborting.\n");
    exit(0);
}

# delete all genes not in this list
# also delete their transcripts, translations and related info
my ($sql, $num);
$support->log_stamped("Deleting internal genes, transcripts and translations...\n");
$sql = qq(
    DELETE QUICK IGNORE
            g,
            gsi,
            gi,
            cgi,
            gn,
            gr,
            gs,
            t,
            tsi,
            ti,
            cti,
            ta,
            tr,
            e,
            tl,
            tla,
            tlsi,
            pf
    FROM
            gene g,
            gene_stable_id gsi,
            gene_info gi,
            current_gene_info cgi,
            gene_name gn,
            transcript t,
            transcript_stable_id tsi,
            transcript_info ti,
            current_transcript_info cti
    LEFT JOIN
            gene_remark gr ON gr.gene_info_id = gi.gene_info_id
    LEFT JOIN
            gene_synonym gs ON gs.gene_info_id = gi.gene_info_id
    LEFT JOIN
            transcript_attrib ta ON ta.transcript_id = t.transcript_id
    LEFT JOIN
            transcript_remark tr ON tr.transcript_info_id = ti.transcript_info_id
    LEFT JOIN
            evidence e ON e.transcript_info_id = ti.transcript_info_id
    LEFT JOIN
            translation tl ON tl.transcript_id = t.transcript_id
    LEFT JOIN
            translation_attrib tla ON tla.translation_id = tl.translation_id
    LEFT JOIN
            translation_stable_id tlsi ON tlsi.translation_id = tl.translation_id
    LEFT JOIN
            protein_feature pf ON pf.translation_id = tl.translation_id
    WHERE   gsi.stable_id NOT IN ($gsi_list)
    AND     g.gene_id = gsi.gene_id
    AND     gsi.stable_id = gi.gene_stable_id
    AND     gsi.stable_id = cgi.gene_stable_id
    AND     gi.gene_info_id = gn.gene_info_id
    AND     t.gene_id = g.gene_id
    AND     t.transcript_id = tsi.transcript_id
    AND     tsi.stable_id = ti.transcript_stable_id
    AND     tsi.stable_id = cti.transcript_stable_id
);
$num = $dbh->do($sql);
$support->log_stamped("Done deleting $num records.\n\n");

# delete exons and supporting features
$support->log_stamped("Deleting exons and supporting features...\n");
$sql = qq(
    DELETE QUICK IGNORE
            e,
            esi,
            et,
            sf
    FROM
            exon e,
            exon_stable_id esi,
            exon_transcript et
    LEFT JOIN
            transcript t ON t.transcript_id = et.transcript_id
    LEFT JOIN
            supporting_feature sf ON sf.exon_id = e.exon_id
    WHERE   t.transcript_id IS NULL
    AND     e.exon_id = esi.exon_id
    AND     e.exon_id = et.exon_id
);
$num = $dbh->do($sql);
$support->log_stamped("Done deleting $num records.\n\n");

# delete xrefs
$support->log_stamped("Deleting xrefs...\n");
$support->log_stamped("Determining which xrefs to delete...\n", 1);
my @xrefs;

# orphan gene xrefs to delete
$sql = qq(
    SELECT
            ox.xref_id,
            ox.object_xref_id
    FROM
            object_xref ox
    LEFT JOIN
            gene g ON g.gene_id = ox.ensembl_id
    WHERE   ox.ensembl_object_type = 'Gene'
    AND     g.gene_id IS NULL
);
my @gene_xrefs = @{ $dbh->selectall_arrayref($sql) || [] };
my $gene_xref_string = join(",", map { $_->[0] } @gene_xrefs);

# since xrefs can be shared between genes, the above list of xrefs might also
# contain entries that are no orphanes, so we have to filter them out
$sql = qq(
    SELECT
            x.xref_id
    FROM
            xref x,
            object_xref ox,
            gene g
    WHERE   g.gene_id = ox.ensembl_id
    AND     ox.ensembl_object_type = 'Gene'
    AND     ox.xref_id = x.xref_id
    AND     x.xref_id IN ($gene_xref_string)
);
my @keep_gene_xrefs = @{ $dbh->selectall_arrayref($sql) || [] };
my %seen_genes;
map { $seen_genes{$_->[0] } = 1 } @keep_gene_xrefs;
foreach my $gene_xref (@gene_xrefs) {
    push(@xrefs, $gene_xref) unless ($seen_genes{$gene_xref->[0]});
}

# orphan transcript xrefs to delete
$sql = qq(
    SELECT
            ox.xref_id,
            ox.object_xref_id
    FROM
            object_xref ox
    LEFT JOIN
            transcript t ON t.transcript_id = ox.ensembl_id
    WHERE   ox.ensembl_object_type = 'Transcript'
    AND     t.transcript_id IS NULL
);
my @transcript_xrefs = @{ $dbh->selectall_arrayref($sql) || [] };
my $transcript_xref_string = join(",", map { $_->[0] } @transcript_xrefs);

# filter (see genes for explanation)
$sql = qq(
    SELECT
            x.xref_id
    FROM
            xref x,
            object_xref ox,
            transcript t
    WHERE   t.transcript_id = ox.ensembl_id
    AND     ox.ensembl_object_type = 'Transcript'
    AND     ox.xref_id = x.xref_id
    AND     x.xref_id IN ($transcript_xref_string)
);
my @keep_transcript_xrefs = @{ $dbh->selectall_arrayref($sql) || [] };
my %seen_transcripts;
map { $seen_transcripts{$_->[0] } = 1 } @keep_transcript_xrefs;
foreach my $transcript_xref (@transcript_xrefs) {
    push(@xrefs, $transcript_xref) unless ($seen_transcripts{$transcript_xref->[0]});
}

# translations xrefs
$sql = qq(
    SELECT
            ox.xref_id,
            ox.object_xref_id
    FROM
            object_xref ox
    LEFT JOIN
            translation tl ON tl.translation_id = ox.ensembl_id
    WHERE   ox.ensembl_object_type = 'Translation'
    AND     tl.translation_id IS NULL
);
my @translation_xrefs = @{ $dbh->selectall_arrayref($sql) || [] };
my $translation_xref_string = join(",", map { $_->[0] } @translation_xrefs);

# filter (see genes for explanation)
$sql = qq(
    SELECT
            x.xref_id
    FROM
            xref x,
            object_xref ox,
            translation tl
    WHERE   tl.translation_id = ox.ensembl_id
    AND     ox.ensembl_object_type = 'Translation'
    AND     ox.xref_id = x.xref_id
    AND     x.xref_id IN ($translation_xref_string)
);
my @keep_translation_xrefs = @{ $dbh->selectall_arrayref($sql) || [] };
my %seen_translations;
map { $seen_translations{$_->[0] } = 1 } @keep_translation_xrefs;
foreach my $translation_xref (@translation_xrefs) {
    push(@xrefs, $translation_xref) unless ($seen_translations{$translation_xref->[0]});
}

my $xref_string = join(",", map { $_->[0] } @xrefs);
my $object_xref_string = join(",", map { $_->[1] } @gene_xrefs, @transcript_xrefs, @translation_xrefs);

# delete from xref
$support->log_stamped("Deleting from xref...\n", 1);
$sql = qq(DELETE FROM xref WHERE xref_id IN ($xref_string));
$num = $dbh->do($sql);
$support->log_stamped("Done deleting $num entries.\n", 1);

# delete from object_xref
$support->log_stamped("Deleting from object_xref...\n", 1);
$sql = qq(DELETE FROM object_xref WHERE object_xref_id IN ($object_xref_string));
$num = $dbh->do($sql);
$support->log_stamped("Done deleting $num entries.\n", 1);

# delete from identity_xref
$support->log_stamped("Deleting from identity_xref...\n", 1);
$sql = qq(DELETE FROM identity_xref WHERE object_xref_id IN ($object_xref_string));
$num = $dbh->do($sql);
$support->log_stamped("Done deleting $num entries.\n", 1);

# delete from external_synonym
$support->log_stamped("Deleting from external_synonym...\n", 1);
$sql = qq(DELETE FROM external_synonym WHERE xref_id IN ($xref_string));
$num = $dbh->do($sql);
$support->log_stamped("Done deleting $num entries.\n", 1);

$support->log_stamped("Done.\n\n");

# optimize tables
$support->log_stamped("Optimizing tables...\n");
my @tables = qw(
    gene
    gene_stable_id
    gene_info
    current_gene_info
    gene_name
    gene_remark
    gene_synonym
    transcript
    transcript_stable_id
    transcript_info
    current_transcript_info
    transcript_attrib
    transcript_remark
    evidence
    translation
    translation_attrib
    translation_stable_id
    protein_feature
    exon
    exon_stable_id
    exon_transcript
    supporting_feature
    xref
    object_xref
    identity_xref
    external_synonym
);
foreach my $table (@tables) {
    $support->log_stamped("$table...\n", 1);
    $dbh->do(qq(OPTIMIZE TABLE $table));
}
$support->log_stamped("Done.\n\n");

# finish logfile
$support->finish_log;

