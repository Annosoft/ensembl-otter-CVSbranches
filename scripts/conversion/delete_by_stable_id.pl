#!/usr/local/bin/perl

=head1 NAME

delete_by_stable_id.pl - filter internal annotation from a Vega db

=head1 SYNOPSIS

delete_by_stable_id.pl [options]

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

    --keep=FILE                         read list of stable IDs to *keep* from
                                        FILE
    --delete=FILE                       read list of stable IDs to *delete*
                                        from FILE
    --find_missing                      print list of genes in infile but not in
                                        database
    --outfile=FILE                      write list of missing genes to FILE
    
=head1 DESCRIPTION

Use this script to delete genes and/or transcripts and all associated
annotation (transcript, translations, xrefs, ...) from a Vega database. You can
provide either a list of gene_stable_ids to keep (--keep=FILE) or to delete
(--delete=FILE), where FILE is a list of stable IDs, one per line. It will
automatically determine whether a stable ID is a gene or a transcript.

The script also checks if any genes/transcripts in the list are not in the
database and optionally prints a list of these stable IDs (use --find_missing
and --outfile=FILE).

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

our $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'keep=s',
    'delete=s',
    'outfile=s',
    'find_missing',
);
$support->allowed_params(
    $support->get_common_params,
    'keep',
    'delete',
    'outfile',
    'find_missing',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# connect to database and get adaptors
my $dba = $support->get_database('ensembl');
our $dbh = $dba->dbc->db_handle;

# find out if we are dealing with an Ensembl or Vega schema
my (%tabs, $schema);
map { $_ =~ s/`//g; $tabs{$_} = 1; } $dbh->tables;
$tabs{'gene_info'} ? $schema = 'vega' : $schema = 'ensembl';

# sanity check: you can only use either a list of gene_stable_id to keep or to
# delete
my ($action, $condition, $infile);
if ($support->param('keep')) {
    $infile = $support->param('keep');
    $action = 'keep';
    $condition = "NOT IN";
} elsif ($support->param('delete')) {
    $infile = $support->param('delete');
    $action = 'delete';
    $condition = "IN";
} else {
    $support->log_error("You must supply either a list of gene_stable_ids to delete or to keep.");
}

# make sure user knows what he's doing
unless ($support->user_proceed("You decided to ".uc($action)." all genes and/or transcripts $condition the list provided. Are you sure you want to proceed?")) {
    exit(0);
}

# read list of stable IDs to keep or delete
my ($gene_stable_ids, $trans_stable_ids) = &read_infile($action, $infile); 

# sanity check: check if all genes in the list are also in the database
&check_missing($gene_stable_ids, $trans_stable_ids);

# nothing else to be done for dry runs
if ($support->param('dry_run')) {
    $support->log("Nothing else to be done for dry run. Aborting.\n");
    exit(0);
}

# delete genes
my $deleted = 0;
$deleted += &delete_genes($gene_stable_ids, $schema, $condition);

# delete transcripts
$deleted += &delete_transcripts($trans_stable_ids, $schema, $condition);

# only try to delete exons and xrefs if you actually deleted genes and/or
# transcripts
if ($deleted) {
    # delete exons
    &delete_exons;

    # delete xrefs
    &delete_xrefs;

    # optimize tables
    &optimize_tables;
}

# finish logfile
$support->finish_log;


### END main ###


=head2 

  Arg[1]      : 
  Example     : 
  Description : read list of stable IDs to keep or delete
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub read_infile {
    my ($action, $infile) = @_;

    $support->log_stamped("Reading stable IDs to ".$action." from file...\n");

    my $in = $support->filehandle('<', $infile);
    my @gene_stable_ids = ();
    my @trans_stable_ids = ();
    while (<$in>) {
        chomp;
        if ($_ =~ /^OTT...G/) {
            push @gene_stable_ids, $_;
        } elsif ($_ =~ /^OTT...T/) {
            push @trans_stable_ids, $_;
        }
    }
    close($in);

    $support->log_stamped("Done reading ".scalar(@gene_stable_ids)." genes, ".scalar(@trans_stable_ids)." transcripts.\n\n");

    return \@gene_stable_ids, \@trans_stable_ids
}

=head2 check_missing

  Arg[1]      : 
  Example     : 
  Description : check if all genes in the list are also in the database
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub check_missing {
    my ($gene_stable_ids, $trans_stable_ids) = @_;

    $support->log("Checking for missing genes and/or transcripts...\n");

    # genes
    my $gsi_string = join("', '", @{ $gene_stable_ids });
    my $sql = qq(
        SELECT  stable_id
        FROM    gene_stable_id
        WHERE   stable_id IN ('$gsi_string')
    );
    my @genes_in_db = map { $_->[0] } @{ $dbh->selectall_arrayref($sql) || [] };
    my $gdiff = scalar(@{ $gene_stable_ids }) - scalar(@genes_in_db);

    # transcripts
    my $tsi_string = join("', '", @{ $trans_stable_ids });
    $sql = qq(
        SELECT  stable_id
        FROM    transcript_stable_id
        WHERE   stable_id IN ('$tsi_string')
    );
    my @trans_in_db = map { $_->[0] } @{ $dbh->selectall_arrayref($sql) || [] };
    my $tdiff = scalar(@{ $trans_stable_ids }) - scalar(@trans_in_db);
    
    if ($gdiff or $tdiff) {
        $support->log_warning("Not all genes and/or transcripts in the input file could be found in the db ($gdiff genes, $tdiff transcripts missing).\n");
        # print list of missing stable IDs
        if ($support->param('find_missing')) {
            $support->log("Printing list of missing stable IDs to file...\n", 1);
            my $out = $support->filehandle('>', $support->param('outfile'));

            # genes
            my %gseen;
            @gseen{@genes_in_db} = (1) x @genes_in_db;
            foreach my $gsi (@{ $gene_stable_ids }) {
                print $out "$gsi\n" unless $gseen{$gsi};
            }

            # transcripts
            my %tseen;
            @tseen{@trans_in_db} = (1) x @trans_in_db;
            foreach my $tsi (@{ $trans_stable_ids }) {
                print $out "$tsi\n" unless $tseen{$tsi};
            }
            
            $support->log("Done.\n", 1);
            
        # suggest to run with --find_missing option
        } else {
            $support->log("Please run the script with the --find_missing option and check to see what's wrong.\n");
        }

        # ask user if he wants to proceed regardless of the potential problem
        unless ($support->user_proceed("\nPotential data inconsistencies (see logfile). Would you like to proceed anyway?")) {
            exit(0);
        }
    }
    $support->log("Done.\n\n");
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 0 if no genes to delete, number of records deleted otherwise
  Exceptions  : 
  Caller      : 

=cut

sub delete_genes {
    my ($gene_stable_ids, $schema, $condition) = @_;

    unless ($gene_stable_ids) {
        $support->log("No genes to delete.\n");
        return(0);
    }

    $support->log_stamped("Deleting genes by stable ID (and their transcripts and translations) ...\n");

    my $gsi_string = join("', '", @{ $gene_stable_ids });
    my $sql;
    if ($schema eq 'vega') {
        # delete statement for Vega schema
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
            WHERE   gsi.stable_id $condition ('$gsi_string')
            AND     g.gene_id = gsi.gene_id
            AND     gsi.stable_id = gi.gene_stable_id
            AND     gsi.stable_id = cgi.gene_stable_id
            AND     gi.gene_info_id = gn.gene_info_id
            AND     t.gene_id = g.gene_id
            AND     t.transcript_id = tsi.transcript_id
            AND     tsi.stable_id = ti.transcript_stable_id
            AND     tsi.stable_id = cti.transcript_stable_id
        );
    } else {
        # delete statement for Ensembl schema
        $sql = qq(
            DELETE QUICK IGNORE
                    g,
                    gsi,
                    t,
                    tsi,
                    ta,
                    tl,
                    tla,
                    tlsi,
                    pf
            FROM
                    gene g,
                    gene_stable_id gsi,
                    transcript t,
                    transcript_stable_id tsi
            LEFT JOIN
                    transcript_attrib ta ON ta.transcript_id = t.transcript_id
            LEFT JOIN
                    translation tl ON tl.transcript_id = t.transcript_id
            LEFT JOIN
                    translation_attrib tla ON tla.translation_id = tl.translation_id
            LEFT JOIN
                    translation_stable_id tlsi ON tlsi.translation_id = tl.translation_id
            LEFT JOIN
                    protein_feature pf ON pf.translation_id = tl.translation_id
            WHERE   gsi.stable_id $condition ('$gsi_string')
            AND     g.gene_id = gsi.gene_id
            AND     t.gene_id = g.gene_id
            AND     t.transcript_id = tsi.transcript_id
        );
    }
    my $num = $dbh->do($sql);
    $support->log_stamped("Done deleting $num records.\n\n");

    return($num);
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 0 if no transcripts to delete, number of records deleted
                otherwise
  Exceptions  : 
  Caller      : 

=cut

sub delete_transcripts {
    my ($trans_stable_ids, $schema, $condition) = @_;

    unless ($trans_stable_ids) {
        $support->log("No transcripts to delete.\n");
        return(0);
    }

    $support->log_stamped("Deleting transcripts by stable ID (and their translations) ...\n");

    my $tsi_string = join("', '", @{ $trans_stable_ids });
    my $sql;
    if ($schema eq 'vega') {
        # delete statement for Vega schema
        $sql = qq(
            DELETE QUICK IGNORE
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
                    transcript t,
                    transcript_stable_id tsi,
                    transcript_info ti,
                    current_transcript_info cti
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
            WHERE   tsi.stable_id $condition ('$tsi_string')
            AND     t.transcript_id = tsi.transcript_id
            AND     tsi.stable_id = ti.transcript_stable_id
            AND     tsi.stable_id = cti.transcript_stable_id
        );
    } else {
        # delete statement for Ensembl schema
        $sql = qq(
            DELETE QUICK IGNORE
                    t,
                    tsi,
                    ta,
                    tl,
                    tla,
                    tlsi,
                    pf
            FROM
                    transcript t,
                    transcript_stable_id tsi
            LEFT JOIN
                    transcript_attrib ta ON ta.transcript_id = t.transcript_id
            LEFT JOIN
                    translation tl ON tl.transcript_id = t.transcript_id
            LEFT JOIN
                    translation_attrib tla ON tla.translation_id = tl.translation_id
            LEFT JOIN
                    translation_stable_id tlsi ON tlsi.translation_id = tl.translation_id
            LEFT JOIN
                    protein_feature pf ON pf.translation_id = tl.translation_id
            WHERE   tsi.stable_id $condition ('$tsi_string')
            AND     t.transcript_id = tsi.transcript_id
        );
    }
    my $num = $dbh->do($sql);
    $support->log_stamped("Done deleting $num records.\n\n");

    # now look for orphan genes and delete them
    $support->log_stamped("Looking for orphan genes...\n");
    $sql = qq(
        SELECT g.gene_id
        FROM gene g
        LEFT JOIN transcript t on g.gene_id = t.gene_id
        WHERE t.gene_id IS NULL;
    );
    my @orphan_genes = map { $_->[0] } @{ $dbh->selectall_arrayref($sql) || [] };

    if (@orphan_genes) {
        $support->log_stamped("Deleting orphan genes...\n", 1);
        
        my $gi_string = join("', '", @orphan_genes);
        if ($schema eq 'vega') {
            # delete statement for Vega schema
            $sql = qq(
                DELETE QUICK IGNORE
                        g,
                        gsi,
                        gi,
                        cgi,
                        gn,
                        gr,
                        gs
                FROM
                        gene g,
                        gene_stable_id gsi,
                        gene_info gi,
                        current_gene_info cgi,
                        gene_name gn
                LEFT JOIN
                        gene_remark gr ON gr.gene_info_id = gi.gene_info_id
                LEFT JOIN
                        gene_synonym gs ON gs.gene_info_id = gi.gene_info_id
                WHERE   g.gene_id $condition ('$gi_string')
                AND     g.gene_id = gsi.gene_id
                AND     gsi.stable_id = gi.gene_stable_id
                AND     gsi.stable_id = cgi.gene_stable_id
                AND     gi.gene_info_id = gn.gene_info_id
            );
        } else {
            # delete statement for Ensembl schema
            $sql = qq(
                DELETE QUICK IGNORE
                        g,
                        gsi
                FROM
                        gene g,
                        gene_stable_id gsi
                WHERE   g.gene_id $condition ('$gi_string')
                AND     g.gene_id = gsi.gene_id
            );
        }
        my $num1 = $dbh->do($sql);
        $num += $num1;
        $support->log_stamped("Done deleting $num1 records.\n", 1);
    } else {
        $support->log("No orphan genes found.\n", 1);
    }
    
    $support->log_stamped("Done.\n\n");

    return($num);
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub delete_exons {
    # delete exons and supporting features
    $support->log_stamped("Deleting exons and supporting features...\n");
    my $sql = qq(
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
    my $num = $dbh->do($sql);
    $support->log_stamped("Done deleting $num records.\n\n");
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub delete_xrefs {
    # delete xrefs
    $support->log_stamped("Deleting xrefs...\n");
    $support->log_stamped("Determining which xrefs to delete...\n", 1);
    
    my ($sql, $num, @xrefs);

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
    my $gene_xref_string = join(",", map { $_->[0] } @gene_xrefs) || 0;

    # since xrefs can be shared between genes, the above list of xrefs might
    # also contain entries that are no orphanes, so we have to filter them out
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
    my $transcript_xref_string = join(",", map { $_->[0] } @transcript_xrefs) || 0;

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
    my $translation_xref_string = join(",", map { $_->[0] } @translation_xrefs) || 0;

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

    my $xref_string = join(",", map { $_->[0] } @xrefs) || 0;
    my $object_xref_string = join(",", map { $_->[1] } @gene_xrefs, @transcript_xrefs, @translation_xrefs) || 0;

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
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub optimize_tables {
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
}

