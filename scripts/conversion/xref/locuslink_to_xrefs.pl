#!/usr/local/bin/perl

=head1 NAME

locuslink_to_xrefs.pl - adds xrefs from LL_tmpl input file

=head1 SYNOPSIS

    locuslink_to_xrefs.pl [options]

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
        --lltmplfile, --lltmpl=FILE         read LL_tmpl from FILE

=head1 DESCRIPTION

This script parses a file downloaded from LocusLink
(ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz) and adds xrefs to
LocusLink and RefSeq_dna. LocusLink identifiers are set as the display name for
the matching genes.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>
Original code by Tim Hubbard <th@sanger.ac.uk>

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
    'lltmplfile|lltmpl=s',
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

# connect to database and get adaptors (caching features on one slice only)
# get an ensembl database for better performance (no otter tables are needed)
my $dba = $support->get_database('ensembl');
$Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::SLICE_FEATURE_CACHE_SIZE = 1;
my $sa = $dba->get_SliceAdaptor();
my $ga = $dba->get_GeneAdaptor();
my $ea = $dba->get_DBEntryAdaptor();
# statement handle for display_xref_id update
my $sth = $dba->dbc->prepare("UPDATE gene SET display_xref_id=? WHERE gene_id=?");

my @gene_stable_ids = $support->param('gene_stable_id');
my %gene_stable_ids = map { $_, 1 } @gene_stable_ids;
my $chr_length = $support->get_chrlength($dba);
my @chr_sorted = $support->sort_chromosomes($chr_length);

# get species name
#my $species = $support->get_species_scientific_name($dba);
my $species = 'Canis familiaris';

# read LL_tmpl file
my ($ll, $lcmap) = &parse_lltmpl($support);

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
    my $num_ll = 0;
    my $num_refseq = 0;
    my $num_case = 0;
    my $num_clone = 0;
    my $num_missing = 0;
    my $gnum = 0;
    foreach my $gene (@$genes) {
        my $gsi = $gene->stable_id;
        my $gid = $gene->dbID;
        my $gene_name = $gene->display_xref->display_id;
        my $lc_gene_name = lc($gene_name);

        # filter to user-specified gene_stable_ids
        if (scalar(@gene_stable_ids)){
            next unless $gene_stable_ids{$gsi};
        }

        $gnum++;
        $support->log("Gene $gene_name ($gid, $gsi)...\n");

        if ($ll->{$gene_name}) {
            # LocusLink xref
            $num_ll++;
            my $dbentry = Bio::EnsEMBL::DBEntry->new(
                    -primary_id => $ll->{$gene_name}->{'LOCUSLINK'},
                    -display_id => $gene_name, 
                    -version    => 1,
                    -release    => 1,
                    -dbname     => 'LocusLink',
            );
            $dbentry->status('KNOWN');
            $gene->add_DBEntry($dbentry);
            unless ($support->param('dry_run')) {
                $ea->store($dbentry, $gid, 'Gene');
                $sth->execute($dbentry->dbID, $gid);
                $support->log("Stored LocusLink xref ".$dbentry->dbID." for gene $gid.\n", 1);
            }

            # RefSeq xref
            if ($ll->{$gene_name}->{'NM'}) {
                $num_refseq++;
                my $dbentry = Bio::EnsEMBL::DBEntry->new(
                        -primary_id => $ll->{$gene_name}->{'NM'},
                        -display_id => $ll->{$gene_name}->{'NM'},
                        -version=>1,
                        -release=>1,
                        -dbname=>"RefSeq_dna",
                );
                $dbentry->status('KNOWNXREF');
                $gene->add_DBEntry($dbentry);
                unless ($support->param('dry_run')) {
                    $ea->store($dbentry, $gid, 'Gene');
                    $support->log("Stored RefSeq_dna xref ".$dbentry->dbID." for gene $gid.\n", 1);
                }
            }
        } elsif ($lcmap->{$lc_gene_name}) {
            # check if case in database might be wrong, by doing lc comparision
            $support->log_warning("Possible case error for $gene_name: ".
                join(',',(@{ $lcmap->{$lc_gene_name} })).".\n", 1);
            $num_case++;
        } elsif ($gene_name =~ /^\w+\.\d+$/ || $gene_name =~ /^\w+\-\w+\.\d+$/) {
            # probably a clone based genename - ok
            $support->log("No LocusLink match for $gene_name (but has clonename based name).\n", 1);
            $num_clone++;
        } else {
            # doesn't look like a clone name, so perhaps mistyped 
            $support->log_warning("No LocusLink match for $gene_name.\n", 1);
            $num_missing++;
        }
    }

    $support->log("\nProcessed $gnum genes (of ".scalar @$genes." on this chromosome).\n");
    $support->log("OK:\n");
    $support->log("LocusLink $num_ll, RefSeq_dna $num_refseq.\n", 1);
    $support->log("WARNINGS:\n");
    $support->log("Genes with possible case mismatch: $num_case.\n", 1);
    $support->log("Genes with apparently clonename based names: $num_clone.\n", 1);
    $support->log("Other genes without match: $num_missing.\n", 1);
    $support->log("Done with chromosome $chr. ".$support->date_and_mem."\n\n");
}

# finish log
$support->log($support->finish_log);

### end main ###

=head2 parse_lltmpl

  Arg[1]      : Bio::EnsEMBL::Utils::ConversionSupport object $support
  Example     : my $ll = &parse_lltmpl($support);
                foreach my $gene (keys %$ll) {
                    my ($nm, $np, $locus_id, $desc) = @{ $ll->{$gene} };
                }
  Description : Parses a LL_tmpl file from LocusLink. File can optionally be
                gzipped.
  Return type : list of hashrefs:
                1.  keys: gene_names
                    values: arrayref of NM, NP, OFFICIAL_SYMBOL, SUMMARY
                2.  keys: lowercase gene_names
                    values: list of gene_names (with case preserved)
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_lltmpl {
    my $support = shift;

    # read LL_tmpl input file
    $support->log("Reading LL_tmpl file... ".$support->date_and_mem."\n");
    $support->check_required_params('lltmplfile');
    my $lltmpl = $support->param('lltmplfile');
    my $fh_expr;
    if($lltmpl =~ /\.gz$/) {
        $fh_expr = "gzip -d -c $lltmpl |";
    } else {
        $fh_expr = "< $lltmpl";
    }
    open(LL, $fh_expr)
        or $support->throw("Couldn't open $lltmpl for reading: $!\n");

    my %ll;
    my %lcmap;
    my $num_total = 1;
    my $num_ok = 0;
    my $num_wrong_org = 0;
    my $num_missing_org = 0;
    my $num_missing_symbol = 0;
    my $locus_id;
    my $flag_found = 1;
    my $flag_org = 0;
    my $nm = '';
    my $np = '';
    my $desc = '';
    my $gene_name = '';

    while(<LL>){
        if (/^SUMMARY: (.*)/) {
            $desc = $1;
        } elsif (/^NM: (.*)/) {
            $nm = $1;
            $nm =~ s/\|.*//;
        } elsif (/^NP: (.*)/) {
            $np = $1;
            $np =~ s/\|.*//;
        } elsif (/^ORGANISM: (.*)/) {
            my $org = $1;
            if ($species eq $org) {
                $flag_org = 1;
            } else {
                $flag_org = 2;
            }
        } elsif (/^OFFICIAL_SYMBOL: (\w+)/) {
            if ($flag_org == 1) {
                $gene_name = $1;
                push @{ $lcmap{lc($gene_name)} }, $gene_name;
                $num_ok++;
            } elsif($flag_org == 2) {
                # wrong organism
                $num_wrong_org++;
            } elsif($flag_org == 0) {
                # missing organism
                $support->log_warning("ORGANISM not found for $locus_id.\n", 1);
                $num_missing_org++;
            }
            $flag_org = 0;
            $flag_found = 1;
        } elsif (/\>\>(\d+)/) {
            if ($gene_name) {
                $ll{$gene_name} = {
                    NM          => $nm,
                    NP          => $np,
                    LOCUSLINK   => $locus_id,
                    DESC        => $desc
                };
                #$support->log("$gene_name\t$nm\t$np\t$locus_id\t$desc\n", 1);
                $nm = '';
                $np = '';
                $desc = '';
            }
            unless ($flag_found) {
                $num_missing_symbol++;
                #$support->log_warning("OFFICIAL_SYMBOL for $locus_id not found.\n", 1);
            }
            $flag_found = 0;
            $locus_id = $1;
            $num_total++;
        }
    }
    close(LL);

    $support->log("Done processing $num_total entries. ".$support->date_and_mem."\n");
    $support->log("OK: $num_ok.\n");
    $support->log("SKIPPED:\n");
    $support->log("Not \'$species\': $num_wrong_org.\n", 1);
    $support->log("No organism label: $num_missing_org.\n", 1);
    $support->log("No official symbol: $num_missing_symbol.\n\n", 1);
    
    return \%ll, \%lcmap;
}

