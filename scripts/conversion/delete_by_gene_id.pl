#!/usr/local/bin/perl

=head1 NAME

delete_by_gene_id.pl - filter internal annotation from a Vega or Ensembl db

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

    --schematype=TYPE                   assume either vega or ensembl schema
    
=head1 DESCRIPTION

Use this script to delete genes and/or transcripts and all associated
annotation (transcript, translations, xrefs, ...) from a Vega or Ensembl database.
The schema type can be provided as a command line argument or determined automaticaly
by the presence or abscence of the gene_info table.

You can provide either a list of gene_stable_ids to keep (--keep=FILE) or to
delete (--delete=FILE), where FILE is a list of stable IDs, one per line. It will
automatically determine whether a stable ID is a gene or a transcript.

The script also checks if any genes/transcripts in the list are not in the
database and optionally prints a list of these stable IDs (use --find_missing
and --outfile=FILE).

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <pm2@sanger.ac.uk>

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
    $support->log_error("You must supply either a list of gene_ids to delete or to keep.");
}

# make sure user knows what he's doing
unless ($support->user_proceed("You decided to ".uc($action)." all genes $condition the list provided. Are you sure you want to proceed?")) {
    exit(0);
}

# read list of gene IDs to keep or delete
my ($gene_ids) = &read_infile($action, $infile); 

# sanity check: check if all genes in the list are also in the database
&check_missing($gene_ids);

# nothing else to be done for dry runs
if ($support->param('dry_run')) {
    $support->log("Nothing else to be done for dry run. Aborting.\n");
    exit(0);
}

# delete genes
my $deleted = 0;
$deleted += &delete_genes($gene_ids, $condition);

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


=head2 read_infile

  Arg[1]      : String $infile - name of file to read gene IDs from
  Example     : my ($gene_ids) = &read_infile('/my/input/file.txt');
  Description : read list of gene IDs to keep or delete from file
  Return type : Arrayref - listref of gene IDs
  Exceptions  : none
  Caller      : internal

=cut

sub read_infile {
    my ($action, $infile) = @_;

    $support->log_stamped("Reading gene IDs to ".$action." from file...\n");

    my $in = $support->filehandle('<', $infile);
    my @gene_ids = ();
    while (<$in>) {
        chomp $_;
		push @gene_ids, $_;
    }
    close($in);
	
    $support->log_stamped("Done reading ".scalar(@gene_ids)." genes.\n\n");

    return \@gene_ids;
}

=head2 check_missing

  Arg[1]      : Arrayref $gene_ids - listref of gene IDs
  Example     : &check_missing($gene_stable_ids);
  Description : Check if all genes in the list are also in the
                database. Warn if this is not the case, and print list of
                missing IDs on request (use --find_missing and --outfile
                options).
  Return type : none
  Exceptions  : none
  Caller      : general

=cut

sub check_missing {
    my ($gene_ids) = @_;

    $support->log("Checking for missing genes...\n");

    # genes
    my $gi_string = join("', '", @{ $gene_ids });
    my $sql = qq(
        SELECT  gene_id
        FROM    gene
        WHERE   gene_id IN ('$gi_string')
    );
    my @genes_in_db = map { $_->[0] } @{ $dbh->selectall_arrayref($sql) || [] };
    my $gdiff = scalar(@{ $gene_ids }) - scalar(@genes_in_db);

    if ($gdiff) {
        $support->log_warning("Not all genes in the input file could be found in the db ($gdiff genes missing).\n");
        # print list of missing stable IDs
        if ($support->param('find_missing')) {
            $support->log("Printing list of missing stable IDs to file...\n", 1);
            my $out = $support->filehandle('>', $support->param('outfile'));

            my %gseen;
            @gseen{@genes_in_db} = (1) x @genes_in_db;
            foreach my $gi (@{ $gene_ids }) {
                print $out "$gi\n" unless $gseen{$gi};
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

=head2 delete_genes

  Arg[1]      : Arrayref $gene_stable_ids - listref of gene stable IDs
  Arg[2]      : String $condition - condition to use in WHERE clause (IN|NOT IN)
  Example     : my $deleted += &delete_genes($gene_stable_ids, 'vega', 'IN');
  Description : Deletes all genes (and associated transcripts and translations)
                in the list provided. 
  Return type : 0 if no genes to delete, number of records deleted otherwise
  Exceptions  : none
  Caller      : internal

=cut

sub delete_genes {
    my ($gene_ids, $condition) = @_;

    unless ($gene_ids) {
        $support->log("No genes to delete.\n");
        return(0);
    }

    $support->log_stamped("Deleting genes by gene ID (and their transcripts and translations) ...\n");

    my $gi_string = join("', '", @{ $gene_ids });
    my $sql= qq(
            DELETE QUICK IGNORE
                    g,
                    gsi,
                    ga,
                    t,
                    tsi,
                    ta,
                    tl,
                    tla,
                    tlsi,
                    pf
            FROM
                    gene g,
                    gene_stable_id gsi
            LEFT JOIN
                    gene_attrib ga ON ga.gene_id = g.gene_id
            LEFT JOIN
                    transcript t on t.gene_id = g.gene_id
            LEFT JOIN
                    transcript_attrib ta ON ta.transcript_id = t.transcript_id
            LEFT JOIN
                    transcript_supporting_feature tsf ON tsf.transcript_id = t.transcript_id
            LEFT JOIN
                    translation tl ON tl.transcript_id = t.transcript_id
            LEFT JOIN
                    translation_attrib tla ON tla.translation_id = tl.translation_id
            LEFT JOIN
                    translation_stable_id tlsi ON tlsi.translation_id = tl.translation_id
            LEFT JOIN
                    protein_feature pf ON pf.translation_id = tl.translation_id
            LEFT JOIN
                    transcript_stable_id tsi on t.transcript_id = tsi.transcript_id
            WHERE   g.gene_id $condition ('$gi_string')
            AND     g.gene_id = gsi.gene_id
        );
    my $num = $dbh->do($sql);
    $support->log_stamped("Done deleting $num records.\n\n");

    return($num);
}


=head2 delete_exons

  Example     : &delete_exons;
  Description : Delete exons (and associated supporting evidence) that don't
                belong to a transcript anymore.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub delete_exons {
    # delete exons and supporting features
	$support->log_stamped("Deleting exon_transcript entries...\n");
	my $sql = qq(
        DELETE QUICK IGNORE
                et
        FROM
                exon_transcript et
        LEFT JOIN
                transcript t on et.transcript_id = t.transcript_id
        WHERE   t.transcript_id is null
    );
	my $num = $dbh->do($sql);
    $support->log_stamped("Done deleting $num records.\n\n");

    $support->log_stamped("Deleting exons and supporting features...\n");
    $sql = qq(
        DELETE QUICK IGNORE
                esi,
                e,
                sf
        FROM
                exon_stable_id esi,
                exon e
        LEFT JOIN
                exon_transcript et ON e.exon_id = et.exon_id
        LEFT JOIN
                supporting_feature sf ON e.exon_id = sf.exon_id
        WHERE   e.exon_id = esi.exon_id
        AND     et.exon_id IS NULL

    );
    $num = $dbh->do($sql);
    $support->log_stamped("Done deleting $num records.\n\n");
}

=head2 delete_xrefs

  Example     : &delete_xrefs;
  Description : Delete xrefs no longer attached to an Ensembl object
  Return type : none
  Exceptions  : none
  Caller      : internal

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

=head2 optimize_tables

  Example     : &optimize_tables;
  Description : Optimises database tables.
  Return type : none
  Exceptions  : none
  Caller      : internal

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
        gene_attrib
        transcript
        transcript_stable_id
        transcript_info
        current_transcript_info
        transcript_attrib
        transcript_remark
        transcript_supporting_feature
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

