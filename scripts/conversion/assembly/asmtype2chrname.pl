#!/usr/local/bin/perl

=head1 NAME

asmtype2chrname.pl - convert assembly type to fake chromosome names

=head1 SYNOPSIS

    asmtype2chrname.pl [options]

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

Otter uses multiple 'type' entries on assembly table. VEGA requires a single
type 'VEGA' instead.

However, its dangerous just to change the type as in cases where there are
multiple regions from a single chromsome in otter, this would cause them to
overlap. One solution is to create a separate chromosome for each type - that's
what this script does.

Pseudocode:
    Get length of region
    create chromosome entry of form num_region
    move assembly of that type to new chromosome
    remove old chromosomes
    move type to VEGA

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
    $SERVERROOT = "$Bin/../../../..";
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
my $dba = $support->get_database('core');
my $dbh = $dba->dbc->db_handle;
my $sql;

# sanity check: make sure you don't have duplicate (ie shared between
# chromosomes/assemblies) contigs in assembly
my $sth = $dbh->prepare('SELECT count(*) FROM assembly');
$sth->execute or $support->throw($sth->errstr);
my ($all) = $sth->fetchrow_array;
$sth->finish;
$sth = $dbh->prepare('SELECT count(distinct(contig_id)) FROM assembly');
$sth->execute or $support->throw($sth->errstr);
my ($unique) = $sth->fetchrow_array;
$sth->finish;
unless ($all == $unique) {
    $support->log("ERROR:\nThere are duplicate (ie shared between chromosomes/assemblies) contigs in the assembly table. You have to fix this before running this script.\nAborting.\n");
    exit(0);
}

# make a backup copy of the chromosome table
$support->log("Backing up chromosome table to chromosome_old...\n");
unless ($support->param('dry_run')) {
    my $b = $dbh->do("CREATE TABLE chromosome_old SELECT * FROM chromosome")
                      or $support->throw($dbh->errstr);
    $support->log("Done backing up $b rows.\n");
}

# get max(chr_length) by assembly type
$support->log("Fetching chromosome.name, assembly.type, chr_length...\n");
my $sth1 = $dbh->prepare(qq(
    SELECT
        c.name          AS chr_name,
        a.type          AS ass_type,
        max(a.chr_end)  AS length
    FROM
        assembly a,
        chromosome c
    WHERE a.chromosome_id = c.chromosome_id
    GROUP by a.type
    ORDER by c.name
));
$sth1->execute or $support->throw($sth1->errstr);
my @rows;
while (my @r = $sth1->fetchrow_array) {
    push @rows, \@r;
}
$support->log("Done fetching ".(scalar @rows)." rows.\n");

# delete old chromosomes
$sql = "DELETE FROM chromosome";
$support->log("Deleting old chromosomes...\n");
$support->log("$sql\n", 1);
unless ($support->param('dry_run')) {
    my $c_rows = $dbh->do($sql) or $support->throw($dbh->errstr);
    $support->log("Done deleting $c_rows rows.\n", 1);
}

# loop over assembly types
$support->log("Looping over chromosomes/assemblies...\n");
my $i = 1;
foreach my $row (@rows) {
    my ($chr, $ass_type, $length) = @$row;
    $support->log("Chr $chr, assembly.type $ass_type...\n", 1);

    # create fake chromosome
    my $sql2 = "INSERT into chromosome values ($i, '$chr-$ass_type', $length)";
    $support->log("$sql2\n", 2);
    unless ($support->param('dry_run')) {
        my $sth2 = $dbh->prepare($sql2);
        $sth2->execute or $support->throw($sth2->errstr);
    }

    # update assembly table with new chromosome names
    my $sql3 = "UPDATE assembly SET chromosome_id = $i WHERE type = '$ass_type'";
    $support->log("$sql3\n", 2);
    unless ($support->param('dry_run')) {
        my $sth3 = $dbh->prepare($sql3);
        $sth3->execute or $support->throw($sth3->errstr);
    }
    
    $i++;
}

# set assembly.type to VEGA for all chromosomes
$sql = qq(UPDATE assembly SET type = 'VEGA');
$support->log("Setting assembly.type to VEGA...\n");
$support->log("$sql\n", 1);
unless ($support->param('dry_run')) {
    my $a_rows = $dbh->do($sql) or $support->throw($dbh->errstr);
    $support->log("Done updating $a_rows rows.\n", 1);
}

# drop backup of chromosome table
$support->log("Drop chromosome backup table...\n");
unless ($support->param('dry_run')) {
    $dbh->do("DROP TABLE chromosome_old") or $support->throw($dbh->errstr);
    $support->log("Done.\n");
}

# finish log
$support->log($support->finish_log);


