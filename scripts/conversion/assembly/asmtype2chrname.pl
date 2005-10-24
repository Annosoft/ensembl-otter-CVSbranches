#!/usr/local/bin/perl

=head1 NAME

asmtype2chrname.pl - convert assembly type to fake chromosome names

=head1 SYNOPSIS

asmtype2chrname.pl [options]

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

=head1 DESCRIPTION

This script does 2 (slightly related) things:

1. Unwanted chromosomes (with all clones, contigs and dangling features) are
deleted from the database.

2. Otter uses multiple 'type' entries on assembly table. VEGA requires a single
type 'VEGA' instead.

However, its dangerous just to change the type as in cases where there are
multiple regions from a single chromosome in otter, this would cause them to
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
$support->allowed_params($support->get_common_params);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# connect to database and get adaptors
my $dba = $support->get_database('vega');
my $dbh = $dba->dbc->db_handle;
my ($sql, $sth);

# delete all chromosomes (clones and dangling features) except those going into
# the release
# show all chromosomes/assemblies
$sql = qq(
    SELECT
            distinct(a.type)    as type,
            c.name              as chr_name,
            a.chromosome_id     as chr_id
    FROM
            assembly a,
            chromosome c
    WHERE   a.chromosome_id = c.chromosome_id
    ORDER BY c.name, a.type;
); 
$sth = $dbh->prepare($sql);
$sth->execute;
my $txt = sprintf "    %-25s%-26s%-15s\n", qw(CHROMOSOME_NAME ASSEMBLY_TYPE CHROMOSOME_ID);
$txt .= "    " . "-"x70 . "\n";
while (my @row = $sth->fetchrow_array) {
    $txt .= sprintf "    %-25s%-26s%-15s\n", $row[1], $row[0], $row[2];
}
$sth->finish;
print "Chromosomes/assemblies currently in the db:\n\n";
print "$txt\n";

# ask user to nominate the ones to delete
print "Which assemblies would you like to DELETE?\n";
print "Enter comma-separated list of assembly_types:\n";
my $ass_to_delete = <>;
chomp $ass_to_delete;
print "\nYou chose to delete these assemblies from the db:\n";
print "$ass_to_delete\n";
exit unless $support->user_proceed("Continue?");

# delete from assembly
if ($ass_to_delete) {
    $ass_to_delete = join(", ", map { $_ =~ s/\s+//g; $dbh->quote($_) } split(",", $ass_to_delete));
    $sql = qq(DELETE FROM assembly WHERE type IN ($ass_to_delete));
    $support->log("Deleting unwanted assemblies...\n");
    $support->log("$sql\n", 1);
    unless ($support->param('dry_run')) {
        my $rows = $dbh->do($sql);
        $support->log("Done deleting $rows rows.\n", 1);
    }
}

# delete orphaned chromosomes
$sql = qq(
    SELECT
            distinct(c.name),
            c.chromosome_id
    FROM
            chromosome c
    LEFT JOIN assembly a
            ON (c.chromosome_id = a.chromosome_id)
    WHERE   a.type is NULL
);
$sth = $dbh->prepare($sql);
$sth->execute;
my @chr_to_del;
while (my @row = $sth->fetchrow_array) {
    push @chr_to_del, $row[1];
}
my $chr_del = join(", ", @chr_to_del); 
if ($chr_del) {
    $sql = qq(DELETE FROM chromosome WHERE chromosome_id IN ($chr_del));
    $support->log("Deleting ".scalar(@chr_to_del)." orphaned chromosomes...\n");
    $support->log("$sql\n", 1);
    unless ($support->param('dry_run')) {
        my $rows = $dbh->do($sql);
        $support->log("Done deleting $rows rows.\n", 1);
    }
}

# delete contigs, clones and dangling features
$sql = qq(
    SELECT
            distinct(ct.clone_id) 
    FROM
            contig ct 
    LEFT JOIN assembly a 
            ON (a.contig_id = ct.contig_id) 
    WHERE a.contig_id IS NULL
);
$sth = $dba->prepare($sql);
$sth->execute;
my $ca = $dba->get_CloneAdaptor;
my @clone_ids;
while (my ($clone_id) = $sth->fetchrow_array) {
    push @clone_ids, $clone_id;
}
$support->log("Deleting ".scalar(@clone_ids)." orphaned clones...\n");

foreach my $clone_id (@clone_ids){
    my $clone = $ca->fetch_by_dbID($clone_id);
    $support->log("Deleting clone ".$clone->embl_id.".".$clone->embl_version."\n", 1);
    $ca->remove($clone) unless ($support->param('dry_run'));
}
$support->log("Done.\n");

# sanity check: make sure you don't have contigs which are shared between
# chromosomes/assemblies
$sth = $dbh->prepare('SELECT count(*) FROM assembly');
$sth->execute;
my ($all) = $sth->fetchrow_array;
$sth->finish;
$sth = $dbh->prepare('SELECT count(distinct(contig_id)) FROM assembly');
$sth->execute;
my ($unique) = $sth->fetchrow_array;
$sth->finish;
unless ($all == $unique) {
    $support->log("ERROR:\nThere are contigs shared between chromosomes/assemblies in the assembly table. You'll have to run assembly_exception.pl to fix this before running this script.\nAborting.\n");
    exit(0);
}

# make a backup copy of the chromosome table
$support->log("Backing up chromosome table to chromosome_old...\n");
unless ($support->param('dry_run')) {
    my $b = $dbh->do("CREATE TABLE chromosome_old SELECT * FROM chromosome");
    $support->log("Done backing up $b rows.\n");
}

# get max(chr_length) by assembly type
$support->log("Fetching chromosome.name, assembly.type, chr_length...\n");
my $sth1 = $dbh->prepare(qq(
    SELECT
        c.name,
        c.chromosome_id,
        a.type,
        max(a.chr_end)
    FROM
        assembly a,
        chromosome c
    WHERE a.chromosome_id = c.chromosome_id
    GROUP by a.type
    ORDER by c.name
));
$sth1->execute;
my @rows;
while (my @r = $sth1->fetchrow_array) {
    push @rows, \@r;
}
$support->log("Done fetching ".(scalar @rows)." rows.\n");
$support->log("These new chromosome names will be used\n");
$support->log("(unless you choose to keep the old ones):\n\n");
$txt = sprintf "    %-20s%-40s\n", "OLD NAME", "NEW NAME";
$txt .= "    " . "-"x70 . "\n";
foreach my $row (@rows) {
    my ($chr, $chr_id, $ass_type) = @$row;
    $txt .= sprintf "    %-20s%-40s\n", $chr, "$chr-$ass_type";
}
$support->log("$txt\n");

# delete old chromosomes
$sql = "DELETE FROM chromosome";
$support->log("Deleting old chromosomes...\n");
$support->log("$sql\n", 1);
unless ($support->param('dry_run')) {
    my $c_rows = $dbh->do($sql);
    $support->log("Done deleting $c_rows rows.\n", 1);
}

# loop over assembly types
$support->log("Looping over chromosomes/assemblies...\n");
foreach my $row (@rows) {
    my ($chr, $chr_id, $ass_type, $length) = @$row;
    $support->log("Chr $chr, assembly.type $ass_type...\n", 1);

    # ask user if he wants to rename chromosome
    my $chr_new = $chr;
    my $input = $support->read_user_input(
        "\nRename $chr to $chr-$ass_type? [y/N/enter other name]"
    );
    if ($input) {
        if ($input eq 'y') {
            $chr_new = "$chr-$ass_type";
        } else {
            $chr_new = $input;
        }
    }

    # create fake chromosome
    my $sql2 = "INSERT into chromosome values ($chr_id, '$chr_new', $length)";
    $support->log("$sql2\n", 2);
    unless ($support->param('dry_run')) {
        my $sth2 = $dbh->prepare($sql2);
        $sth2->execute;
    }

    # update assembly table with new chromosome names
    my $sql3 = "UPDATE assembly SET chromosome_id = $chr_id WHERE type = '$ass_type'";
    $support->log("$sql3\n", 2);
    unless ($support->param('dry_run')) {
        my $sth3 = $dbh->prepare($sql3);
        $sth3->execute;
    }
}

# set assembly.type to VEGA for all chromosomes
$sql = qq(UPDATE assembly SET type = 'VEGA');
$support->log("Setting assembly.type to VEGA...\n");
$support->log("$sql\n", 1);
unless ($support->param('dry_run')) {
    my $a_rows = $dbh->do($sql);
    $support->log("Done updating $a_rows rows.\n", 1);
}

# add assembly.default value to meta table
$sql = qq(INSERT INTO meta (meta_key, meta_value) VALUES ('assembly.default', 'VEGA'));
$support->log("Adding assembly.default to meta table...\n");
$support->log("$sql\n", 1);
unless ($support->param('dry_run')) {
    my $a_rows = $dbh->do($sql);
    $support->log("Done.\n", 1);
}

# drop backup of chromosome table
$support->log("Drop chromosome backup table...\n");
unless ($support->param('dry_run')) {
    $dbh->do("DROP TABLE chromosome_old");
    $support->log("Done.\n");
}

# finish log
$support->finish_log;

