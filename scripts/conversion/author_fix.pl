#!/usr/local/bin/perl

=head1 NAME

author_fix.pl - sets the right author for genes/transcripts

=head1 SYNOPSIS

    author_fix.pl [options]

    General options:
        --dbname, db_name=NAME              use database NAME
        --host, --dbhost, --db_host=HOST    use database host HOST
        --port, --dbport, --db_port=PORT    use database port PORT
        --user, --dbuser, --db_user=USER    use database username USER
        --pass, --dbpass, --db_pass=PASS    use database passwort PASS
        --driver, --dbdriver, --db_driver=DRIVER    use database driver DRIVER
        --conffile, --conf=FILE             read parameters from FILE
        --logfile, --log=FILE               log to FILE (default: *STDOUT)
        -i, --interactive                   run script interactively
                                            (default: true)
        -n, --dry_run, --dry                don't write results to database
        -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

In Otter, the author of genes/transcripts corresponds to the individual
annotator how created the object. On the Vega website, we display the
annotation group who annotated a particular chromsosome/dataset as the author.
Furthermore, the author is used as a key to group genesets in contigview.

This script set the correct author for genes/transcripts. It uses a lookup hash
to determine which group to use for a chromsome.

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

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->log_filehandle('>>');
$support->log($support->init_log);

# connect to database and get adaptors
my $dba = $support->get_database('otter');
my $dbh = $dba->dbc->db_handle;

# author lookup hash (uses taxonomy_id for species)
my $author_def = {
    # Homo_sapiens
    '9609'  => {
        'default'   => [ 'Havana', 'vega@sanger.ac.uk' ],
        'other'     => {
            '7'         => [ 'Washu', 'jspieth@watson.wust' ],
            '14'        => [ 'Genoscope', 'ralph@genoscope.cns.fr' ],
            '22'        => [ 'Sanger', 'chr22@sanger.ac.uk' ],
        },
    },
    # Mus_musculus
    '10090'  => { 
        'default'    => [ 'Havana', 'vega@sanger.ac.uk' ],
    },
    # Danio_rerio
    '7955'   => {
        'default'    => [ 'zfish', 'zfish-help@sanger.ac.uk' ],
    },
    # Canis_familiaris
    '9615'   => {
        'default'    => [ 'zfish', 'zfish-help@sanger.ac.uk' ],
    },
};

# determine species from database
my $tid = $support->get_taxonomy_id($dba);
my $chromosomes = $author_def->{$tid};
$support->throw("Unknown taxonomy_id. Please update the definition hash in this script.") unless $chromosomes;

# ask user to confirm author lookup hash
print "These author settings will be used:\n\n";
printf "    %-20s%-12s%-34s\n", qw(CHROMOSOME AUTHOR EMAIL);
print "    " . "-"x70 . "\n";
printf "    %-20s%-12s%-34s\n", 'default', @{ $chromosomes->{'default'} };
foreach my $chr (sort keys %{ $chromosomes->{'other'} } ) {
    printf "    %-20s%-12s%-34s\n", $chr, @{ $chromosomes->{'other'}->{$chr} };
}
print "\n";
exit unless $support->user_proceed("Continue?");

if ($support->param('dry_run')) {
    $support->log("There is nothing else to do for dry_run. Aborting.\n");
    exit(0);
}

# delete old authors, insert new authors into database
$support->log("Deleting old authors, inserting new authors into db...\n");
$dbh->do('DELETE FROM author');
my ($author, $email) = @{ $chromosomes->{'default'} };
$dbh->do("INSERT INTO author VALUES (1000, '$email', '$author')");
my $i = 1001;
foreach my $chr (keys %{ $chromosomes->{'other'} }) {
    my ($author, $email) = @{ $chromosomes->{'other'}->{$chr} };
    $dbh->do("INSERT INTO author VALUES ($i, '$email', '$author')");
    $i++;
}
$support->log("Done.\n");

# if we have authors other than the default, go the complicated way
if (scalar(keys %{ $chromosomes->{'other'} })) {
    # Create temporary tables mapping gene_info_id's and transcript_info_id's
    # to chromosome
    $support->log("Creating temporary table with gene_info_id -> chromosome mappings...\n");
    $dbh->do(qq(
        CREATE TABLE chr_gene_info_temp
        SELECT DISTINCT(gi.gene_info_id), sr.name
        FROM
            gene g,
            gene_stable_id gsi,
            current_gene_info cgi,
            gene_info gi,
            seq_region sr
        WHERE   g.gene_id = gsi.gene_id
            AND gsi.stable_id = cgi.gene_stable_id
            AND cgi.gene_info_id = gi.gene_info_id
            AND g.seq_region_id = sr.seq_region_id
    ));
    $support->log("Done.\n");

    $support->log("Creating temporary table with transcript_info_id -> chromosome mappings...\n");
    $dbh->do(qq(
        CREATE TABLE chr_transcript_info_temp
        SELECT DISTINCT(ti.transcript_info_id), sr.name
        FROM
            transcript t,
            transcript_stable_id tsi,
            current_transcript_info cti,
            transcript_info ti,
            seq_region sr
        WHERE   t.transcript_id = tsi.transcript_id
            AND tsi.stable_id = cti.transcript_stable_id
            AND cti.transcript_info_id = ti.transcript_info_id
            AND t.seq_region_id = sr.seq_region_id
    ));
    $support->log("Done.\n");
}

# set all author to default for all genes/transcripts
$support->log("Setting author to default for all genes/transcripts...\n");
$dbh->do(qq(UPDATE gene_info SET author_id = 1000));
$dbh->do(qq(UPDATE transcript_info SET author_id = 1000));
$support->log("Done.\n");

# set correct author_id for chromosomes with non-default author
$support->log("Setting author for chromosomes with non-default author...\n");
foreach my $chr (keys %{ $chromosomes->{'other'} }) {
    my ($author, $email) = @{ $chromosomes->{'other'}->{$chr} };
    $dbh->do(qq(
        UPDATE gene_info gi, chr_gene_info_temp cg, author au
        SET gi.author_id = au.author_id
        WHERE   cg.gene_info_id = gi.gene_info_id
            AND cg.name = '$chr'
            AND au.author_name = '$author'
    ));
    $dbh->do(qq(
        UPDATE transcript_info ti, chr_transcript_info_temp ct, author au
        SET ti.author_id = au.author_id
        WHERE   ct.transcript_info_id = ti.transcript_info_id
            AND ct.name = '$chr'
            AND au.author_name = '$author'
    ));
}
$support->log("Done.\n");

# drop temporary tables
$support->log("Dropping temporary tables...\n");
$dbh->do('DROP TABLE IF EXISTS chr_gene_info_temp');
$dbh->do('DROP TABLE IF EXISTS chr_transcript_info_temp');
$support->log("Done.\n");

# finish logfile
$support->log($support->finish_log);

