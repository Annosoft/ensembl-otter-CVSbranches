#!/usr/local/bin/perl

=head1 NAME

accession_to_support.pl - script to add supporting evidence to a Vega database

=head1 SYNOPSIS

accession_to_support.pl [options]

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
    --chromosomes, --chr=LIST           only process LIST chromosomes
    --gene_stable_id, --gsi=LIST|FILE   only process LIST gene_stable_ids
                                        (or read list from FILE)

=head1 DESCRIPTION

This script adds the supporting evidence for Vega. It does so by comparing
accessions between annotated evidence and similarity features from the protein
pipeline run. If a match is found, it is added to the supporting_feature
table.

Pseudocode:

    foreach gene
        get all similarity features, store in datastructure
        foreach transcript
            get all annotated evidence
            foreach evidence
                foreach similarity feature
                    accession matches?
                        foreach exon
                            similarity feature overlaps exon?
                                store supporting evidence

There are occasions where no match for annotated evidence can be found.
Possible reasons for this are: spelling mistake by annotator; feature not found
by protein pipeline run (e.g. removed from external database, renamed)

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
);
$support->allowed_params(
    $support->get_common_params,
    'chromosomes',
    'gene_stable_id',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');
$support->list_or_file('gene_stable_id');

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# connect to database and get adaptors (caching features on one slice only)
my $dba = $support->get_database('otter');
$Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::SLICE_FEATURE_CACHE_SIZE = 1;
my $sa = $dba->get_SliceAdaptor();
my $ga = $dba->get_GeneAdaptor();
my $aa = $dba->get_AnalysisAdaptor();
# statement handle for storing supporting evidence
my $sql = "insert into supporting_feature
           (exon_id, feature_id, feature_type)
           values(?, ?, ?)";
my $sth = $dba->dbc->prepare($sql);

my @gene_stable_ids = $support->param('gene_stable_id');
my %gene_stable_ids = map { $_, 1 } @gene_stable_ids;
my $chr_length = $support->get_chrlength($dba);
my @chr_sorted = $support->sort_chromosomes($chr_length);
my %analysis = map { $_->logic_name => $_ } @{ $aa->fetch_all };
my %ftype = (
    'Bio::EnsEMBL::DnaDnaAlignFeature' => 'dna_align_feature',
    'Bio::EnsEMBL::DnaPepAlignFeature' => 'protein_align_feature',
);

# loop over chromosomes
$support->log("Looping over chromosomes: @chr_sorted\n\n");
foreach my $chr (@chr_sorted) {
    $support->log("> Chromosome $chr (".$chr_length->{$chr}
               ."bp). ".$support->date_and_mem."\n\n");
    
    # fetch genes from db
    $support->log("Fetching genes...\n");
    my $slice = $sa->fetch_by_region('chromosome', $chr);
    my $genes = $ga->fetch_by_Slice($slice);
    $support->log("Done fetching ".scalar @$genes." genes. " .
                   $support->date_and_mem."\n\n");

    # loop over genes
    my %stats = map { $_ => 0 } qw(genes transcripts exons genes_without_support transcripts_without_support);
    my @transcripts_without_support;
    foreach my $gene (@$genes) {
        my $gsi = $gene->stable_id;
        my $gid = $gene->dbID;
        my $gene_name = $gene->display_xref->display_id;

        # filter to user-specified gene_stable_ids
        if (scalar(@gene_stable_ids)){
            next unless $gene_stable_ids{$gsi};
        }

        # adjust gene's slice to cover gene +/- 1000 bp
        my $gene_slice = $sa->fetch_by_region('chromosome', $chr, $gene->start - 1000, $gene->end + 1000);
        $gene = $gene->transfer($gene_slice);
        unless ($gene) {
            $support->log_warning("Gene $gene_name ($gid, $gsi) doesn't transfer to padded gene_slice.\n");
            next;
        }
        
        $stats{'genes'}++;
        my %se_hash = ();
        my $gene_has_support = 0;
        $support->log_verbose("Gene $gene_name ($gid, $gsi) on slice ".
                        $gene->slice->name."... ".
                        $support->date_and_mem."\n");

        # fetch similarity features from db and store required information in
        # lightweight datastructure (name => [ start, end, dbID, type ])
        $support->log_verbose("Fetching similarity features... ".
                        $support->date_and_mem."\n", 1);
        my $similarity = $gene_slice->get_all_SimilarityFeatures;
        my $sf = {};
        foreach my $f (@$similarity) {
            (my $hitname = $f->hseqname) =~ s/\.[0-9]*$//;
            push @{ $sf->{$hitname} },
                 [ $f->start, $f->end, $f->dbID, $ftype{ref($f)} ];
        }
        $support->log_verbose("Done fetching ".(scalar @$similarity).
                        " features.".$support->date_and_mem."\n", 1);

        # loop over transcripts
        foreach my $trans (@{ $gene->get_all_Transcripts }) {
            my $transcript_has_support = 0;
            $stats{'transcripts'}++;
            $support->log_verbose("Transcript ".$trans->stable_id."...\n", 1);

            # loop over evidence added by annotators for this transcript
            my @evidence = $trans->transcript_info->evidence;
            my @exons = @{ $trans->get_all_Exons };
            $stats{'exons'} += scalar(@exons);
            foreach my $evi (@evidence) {
                my $acc = $evi->name;
                $acc =~ s/.*://;
                $acc =~ s/\.[0-9]*$//;
                my $ana = $analysis{$evi->type . "_evidence"};
                $support->log_verbose("Evidence $acc...\n", 2);
                # loop over similarity features on the slice, compare name with
                # evidence
                my $match = 0;
                foreach my $hitname (keys %$sf) {
                    if ($hitname eq $acc) {
                        foreach my $hit (@{ $sf->{$hitname} }) {
                            # loop over exons and look for overlapping
                            # similarity features
                            foreach my $exon (@exons) {
                                if ($exon->end >= $hit->[0] && $exon->start <= $hit->[1]) {
                                    $support->log_verbose("Matches similarity feature with dbID ".$hit->[2].".\n", 3);
                                    # store unique evidence identifier in hash
                                    $se_hash{$exon->dbID.":".$hit->[2].":".$hit->[3]} = 1;
                                    
                                    $match = 1;
                                    $gene_has_support++;
                                    $transcript_has_support++;
                                }
                            }
                        }
                    }
                }
                $support->log_warning("No matching similarity feature found for $acc.\n", 3) unless ($match);
            }
            unless ($transcript_has_support) {
                $stats{'transcripts_without_support'}++;
                push @transcripts_without_support, $trans->stable_id." on gene ".$gsi;
            }
        }
        $stats{'genes_without_support'}++ unless ($gene_has_support);

        $support->log_verbose("Found $gene_has_support matches (".
                       scalar(keys %se_hash)." unique).\n", 1);

        # store supporting evidence in db
        if ($gene_has_support and !$support->param('dry_run')) {
            $support->log_verbose("Storing supporting evidence... ".
                           $support->date_and_mem."\n", 1);
            foreach my $se (keys %se_hash) {
                eval {
                    $sth->execute(split(":", $se));
                };
                if ($@) {
                    $support->log_warning("$gsi: $@\n", 1);
                }
            }
            $support->log_verbose("Done storing evidence. ".
                           $support->date_and_mem."\n", 1);
        }
    }
    $support->log("\nProcessed $stats{genes} genes (of ".scalar @$genes." on chromosome $chr), $stats{transcripts} transcripts, $stats{exons} exons.\n");
    $support->log("WARNINGS:\n");
    if ($stats{'genes_without_support'}) {
        $support->log("No supporting evidence for any transcripts on $stats{genes_without_support} genes.\n", 1);
        $support->log("No supporting evidence for $stats{transcripts_without_support} transcripts.\n", 1);
        $support->log("Transcripts without supporting evidence:\n", 1);
        foreach (@transcripts_without_support) {
            $support->log("$_\n", 2);
        }
    } else {
        $support->log("None.\n");
    }
    $support->log("Done with chromosome $chr. ".$support->date_and_mem."\n\n");
}

# finish log
$support->finish_log;

