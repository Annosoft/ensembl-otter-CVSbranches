#!/usr/local/bin/perl

=head1 NAME

hugo_to_xrefs.pl - adds xrefs from HUGO input file

=head1 SYNOPSIS

    hugo_to_xrefs.pl [options]

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

    Specific options:
        --chromosomes, --chr=LIST           only process LIST chromosomes
        --gene_stable_id, --gsi=LIST|FILE   only process LIST gene_stable_ids
                                            (or read list from FILE)
        --nomeidfile, --nome=FILE           read HUGO nomenclature from FILE

=head1 DESCRIPTION

This script parses a file downloaded from HUGO
(http://www.gene.ucl.ac.uk/nomenclature/) and adds xrefs to HUGO, LocusLink,
Swissprot, RefSeq, OMIM and others. HUGO names are set as the display name for
the matching genes.

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
$support->parse_extra_options(
    'chromosomes|chr=s@',
    'gene_stable_id|gsi=s@',
    'nomeidfile|nome=s',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');
$support->list_or_file('gene_stable_id');

# ask user to confirm parameters to proceed
$support->confirm_params;

# make sure add_vega_xrefs.pl has been run
print "This script must run after add_vega_xrefs.pl. Have you run it?\n";
$support->user_confirm;

# get log filehandle and print heading and parameters to logfile
$support->log_filehandle('>>');
$support->log($support->init_log);

# read HUGO data from file
$support->log("Reading HUGO nomenclature file...\n");
$support->check_required_params('nomeidfile');
open (NOM, '<', $support->param('nomeidfile')) or die
    "Couldn't open ".$support->param('nomeidfile')."for reading: $!\n";
my $line = <NOM>;
chomp $line;
my @fieldnames = split /\t/, $line;
my $num_fields = scalar(@fieldnames);
my %hugohash;
while (<NOM>) {
    chomp;
    my @fields = split /\t/, $_, -1;
    if (scalar(@fields) != $num_fields) {
        $support->log("Wrong number of fields (got " . scalar(@fields) .
                      ", should be $num_fields).\n");
        $support->log("Entry:\n$_\n");
        $support->log("Aborting.\n");
        die("Inconsistent field numbers in HUGO file. Aborting.\n");
    }
    if (defined($hugohash{$fields[1]})) {
        $support->log("Duplicate lines for " . $fields[1] . ". Aborting.\n");
        die("Duplicate lines in HUGO file. Aborting.\n");
    }
    $hugohash{$fields[1]} = \@fields;
}
$support->log("Done reading ".(scalar keys %hugohash)." entries. ".$support->date_and_mem."\n");

my %convhash = (
    'MIM'           => 'OMIM',
    'Ref Seq'       => 'RefSeq',
    'Locus Link'    => 'LocusLink',
    'SWISSPROT'     => 'SWISSPROT',
);

# connect to database and get adaptors (caching features on one slice only)
# get an ensembl database for better performance (no otter tables are needed)
my $dba = $support->get_database('ensembl');
$Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::SLICE_FEATURE_CACHE_SIZE = 1;
my $sa = $dba->get_SliceAdaptor();
my $ga = $dba->get_GeneAdaptor();
my $ea = $dba->get_DBEntryAdaptor();
# statement handle for display_xref_id update
my $sth = $dba->dbc->prepare("update gene set display_xref_id=?  where gene_id=?");

my @gene_stable_ids = $support->param('gene_stable_id');
my %gene_stable_ids = map { $_, 1 } @gene_stable_ids;
my $chr_length = $support->get_chrlength($dba);
my @chr_sorted = $support->sort_chromosomes($chr_length);

# loop over chromosomes
$support->log("Looping over chromosomes: @chr_sorted\n\n");
foreach my $chr (@chr_sorted) {
    $support->log("> Chromosome $chr (".$chr_length->{$chr}
               ."bp). ".$support->date_and_mem."\n\n");
    
    # fetch genes from db
    $support->log("Fetching genes...\n");
    my $slice = $sa->fetch_by_region('chromosome', $chr);
    my $genes = $ga->fetch_all_by_Slice($slice);
    $support->log("Done fetching ".scalar @$genes." genes. " .
                   $support->date_and_mem."\n\n");

    # loop over genes
    my $n_hugo = 0;
    my %n_other = map { $convhash{$_}, 0 } keys %convhash;
    my $gnum = 0;
    foreach my $gene (@$genes) {
        my $gsi = $gene->stable_id;
        my $gid = $gene->dbID;
        my $gene_name = $gene->display_xref->display_id;

        # filter to user-specified gene_stable_ids
        if (scalar(@gene_stable_ids)){
            next unless $gene_stable_ids{$gsi};
        }

        $gnum++;
        $support->log("Gene $gene_name ($gid, $gsi)...\n");

        # Human hugo symbols are meant to be upper case apart from orfs.
        # There's one which isn't (IL27w).
        my $uc_gene_name;
        if ($gene_name =~ /C.*orf[0-9]*/) {
            $uc_gene_name = $gene_name;
        } else {
            $uc_gene_name = uc $gene_name;
        }

        if (defined($hugohash{$uc_gene_name})) {
            $n_hugo++;
            my $dbentry = Bio::EnsEMBL::DBEntry->new(
                    -primary_id => $hugohash{$uc_gene_name}->[0],
                    -display_id => $gene_name, 
                    -version    => 1,
                    -release    => 1,
                    -dbname     => 'HUGO',
            );
            $dbentry->status('KNOWN');
            $gene->add_DBEntry($dbentry);
            unless ($support->param('dry_run')) {
                $ea->store($dbentry, $gid, 'Gene');
                $sth->execute($dbentry->dbID, $gid);
                $support->log("Stored HUGO xref ".$dbentry->dbID." for gene $gid.\n", 1);
            }

            # store other xrefs for this gene
            for (my $i = 4; $i < 13; $i++) {
                my $xid = $hugohash{$uc_gene_name}->[$i];
                if (exists($convhash{$fieldnames[$i]}) and $xid ne "") {
                    $n_other{$fieldnames[$i]}++;
                    my $dbentry = Bio::EnsEMBL::DBEntry->new(
                            -primary_id => $xid,
                            -display_id => $xid, 
                            -version    => 1,
                            -release    => 1,
                            -dbname     => $convhash{$fieldnames[$i]},
                    );
                    if ($fieldnames[$i] eq "SWISSPROT") {
                        $dbentry->status('XREF');
                    } else {
                        $dbentry->status('KNOWNXREF');
                    }
                    $gene->add_DBEntry($dbentry);
                    unless ($support->param('dry_run')) {
                        $ea->store($dbentry, $gid, 'Gene');
                        $sth->execute($dbentry->dbID, $gid);
                        $support->log("Stored ".$fieldnames[$i]." xref ".$dbentry->dbID." ($xid) for gene $gid.\n", 1);
                    }
                }
            }
        } else {
            $support->log("No HUGO match for gene $gene_name.\n", 1);
        }
    }
    $support->log("\nProcessed $gnum genes (of ".scalar @$genes." on this chromosome).\n");
    $support->log("Genes with no HUGO match: ".($gnum - $n_hugo).".\n");
    my $matches = join ", ", map { "$convhash{$_} $n_other{$_}" } sort keys %convhash;
    $support->log("Matches: HUGO $n_hugo, $matches.\n");
    $support->log("Done with chromosome $chr. ".$support->date_and_mem."\n\n");
}

# finish log
$support->log($support->finish_log);

