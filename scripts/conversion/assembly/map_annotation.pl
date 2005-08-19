#!/usr/local/bin/perl

=head1 NAME

map_annotation.pl - map features from one assembly onto
another

=head1 SYNOPSIS

map_annotation.pl [options]

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
    --evegadbname=NAME                  use ensembl-vega (target) database NAME
    --evegahost=HOST                    use ensembl-vega (target) database host
                                        HOST
    --evegaport=PORT                    use ensembl-vega (target) database port
                                        PORT
    --evegauser=USER                    use ensembl-vega (target) database
                                        username USER
    --evegapass=PASS                    use ensembl-vega (target) database
                                        passwort PASS

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
    $SERVERROOT = "$Bin/../../../..";
    unshift(@INC, "./modules");
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;
use InterimTranscript;
use InterimExon;
use StatMsg;
use Deletion;
use Transcript;
use Gene;
use StatLogger;
use StatMsg;
use Utils qw(print_exon print_coords print_translation);

$| = 1;

our $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'evegahost=s',
    'evegaport=s',
    'evegauser=s',
    'evegapass=s',
    'evegadbname=s',
    'statlog=s',
    'chromosomes|chr=s@',
    'prune=s',
);
$support->allowed_params(
    $support->get_common_params,
    'evegahost',
    'evegaport',
    'evegauser',
    'evegapass',
    'evegadbname',
    'statlog',
    'chromosomes',
    'prune',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# check required params
$support->check_required_params(qw(statlog));

# connect to database and get adaptors
my $V_dba = $support->get_database('core');
my $V_dbh = $V_dba->dbc->db_handle;
my $V_sa = $V_dba->get_SliceAdaptor;
my $V_ga = $V_dba->get_GeneAdaptor;
my $V_pfa = $V_dba->get_ProteinFeatureAdaptor;
my $E_dba = $support->get_database('evega', 'evega');
my $E_dbh = $E_dba->dbc->db_handle;
my $E_sa = $E_dba->get_SliceAdaptor;
my $E_ga = $E_dba->get_GeneAdaptor;
my $E_pfa = $E_dba->get_ProteinFeatureAdaptor;

StatMsg::set_logger(StatLogger->new(
    $support->param('logpath').'/'.$support->param('statlog')));

my $cs_adaptor = $E_dba->get_CoordSystemAdaptor;
my $asmap_adaptor = $E_dba->get_AssemblyMapperAdaptor;

my $E_cs = $cs_adaptor->fetch_by_name('chromosome',
    $support->param('ensemblassembly'));
my $V_cs = $cs_adaptor->fetch_by_name('chromosome',
    $support->param('assembly'));

my $mapper = $asmap_adaptor->fetch_by_CoordSystems($E_cs, $V_cs);
$mapper->max_pair_count( 6_000_000 );
$mapper->register_all;

if ($support->param('prune')) {
    $support->log("Deleting db entries from previous runs of this script...\n");
    $E_dbh->do(qq(DELETE FROM analysis));
    $E_dbh->do(qq(DELETE FROM dna_align_feature));
    $E_dbh->do(qq(DELETE FROM exon));
    $E_dbh->do(qq(DELETE FROM exon_stable_id));
    $E_dbh->do(qq(DELETE FROM exon_transcript));
    $E_dbh->do(qq(DELETE FROM gene));
    $E_dbh->do(qq(DELETE FROM gene_stable_id));
    $E_dbh->do(qq(DELETE FROM interpro));
    $E_dbh->do(qq(DELETE FROM object_xref));
    $E_dbh->do(qq(DELETE FROM protein_align_feature));
    $E_dbh->do(qq(DELETE FROM protein_feature));
    $E_dbh->do(qq(DELETE FROM supporting_feature));
    $E_dbh->do(qq(DELETE FROM transcript));
    $E_dbh->do(qq(DELETE FROM transcript_stable_id));
    $E_dbh->do(qq(DELETE FROM translation));
    $E_dbh->do(qq(DELETE FROM translation_stable_id));
    $E_dbh->do(qq(
        DELETE x
        FROM xref x, external_db ed
        WHERE x.external_db_id = ed.external_db_id
        AND ed.db_name NOT IN ('Vega_gene','Vega_transcript','Vega_translation')
    ));
    $support->log("Done.\n");
}

# loop over chromosomes
$support->log("Looping over chromosomes...\n");
my $V_chrlength = $support->get_chrlength($E_dba, $support->param('assembly'));
my $E_chrlength = $support->get_chrlength($E_dba, $support->param('ensemblassembly'));
foreach my $chr ($support->sort_chromosomes($V_chrlength)) {
    $support->log_stamped("Chromosome $chr...\n", 1);
    
    # skip non-ensembl chromosomes (e.g. MHC haplotypes)
    unless ($E_chrlength->{$chr}) {
        $support->log("Chromosome not in Ensembl. Skipping.\n", 1);
        next;
    }

    # fetch chromosome slices
    my $V_slice = $V_sa->fetch_by_region('chromosome', $chr, undef, undef,
        undef, $support->param('assembly'));
    my $E_slice = $E_sa->fetch_by_region('chromosome', $chr, undef, undef,
        undef, $support->param('ensemblassembly'));

    $support->log("Looping over genes...\n", 1);
    my $genes = $V_ga->fetch_all_by_Slice($V_slice);
    foreach my $gene (@{ $genes }) {
        $support->log("Gene ".$gene->stable_id."\n", 2);

        my $transcripts = $gene->get_all_Transcripts;
        my (@finished, %all_protein_features);
        foreach my $transcript (@{ $transcripts }) {
            my $interim_transcript = transfer_transcript($transcript, $mapper,
                $V_cs, $V_pfa, $E_slice);
            my ($finished_transcripts, $protein_features) =
                create_transcripts($interim_transcript, $E_sa);

            # set the translation stable identifier on the finished transcripts
            foreach my $tr (@{ $finished_transcripts }) {
                if ($tr->translation && $transcript->translation) {
                    $tr->translation->stable_id($transcript->translation->stable_id);
                    $tr->translation->version($transcript->translation->version);
                }
            }

            # This method call is optional, just for statistics gathering
            # and extra sanity checking / debug output:
            #gen_transcript_stats($finished_transcripts);

            push @finished, @$finished_transcripts;
            map { $all_protein_features{$_} = $protein_features->{$_} }
                keys %{ $protein_features || {} };
        }

        unless ($support->param('dry_run')) {
            Gene::store_gene($support, $E_slice, $E_ga, $E_pfa, $gene,
                \@finished, \%all_protein_features);
        }
    }
    $support->log("Done with chromosome $chr.\n", 1);
}
$support->log("Done.\n");

# finish logfile
$support->finish_log;


### END main


#
# method which does some sanity checking of generated transcripts,
# gathers some stats about the completed transcripts and prints the
# created translations to STDERR
#
sub gen_transcript_stats {
  my $finished_transcripts = shift;

  my $transcript_count = @$finished_transcripts;
  my $translation_count = 0;
  my $stop_codons_count = 0;

  if($transcript_count > 1) {
    StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::SPLIT);
  }
  elsif($transcript_count== 0) {
    StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::NO_SEQUENCE_LEFT);
  }

  foreach my $ftrans (@$finished_transcripts) {
    if($ftrans->translation) {
      $translation_count++;
      my $pep = $ftrans->translate->seq;

      print STDERR "\n\n$pep\n\n";

      if($pep =~ /\*/) {
        $stop_codons_count++;
      }

      # sanity check, if translation is defined we expect a peptide
      if(!$pep) {
        print_translation($support, $ftrans->translation);
        $support->log_error("Unexpected Translation but no peptide.\n");
      }
    }
  }

  # If there were stop codons in one of the split transcripts
  # report it. Report it as 'entire' if all split transcripts had
  # stops.
  if($stop_codons_count) {
    my $code = StatMsg::TRANSCRIPT | StatMsg::DOESNT_TRANSLATE;
    if($stop_codons_count == $translation_count) {
      $code |= StatMsg::ENTIRE;
    } else {
      $code |= StatMsg::PARTIAL;
    }
    StatMsg->new($code);
  }

  if(!$translation_count) {
    StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::NO_CDS_LEFT);
  }

  if($translation_count) {
    if($stop_codons_count) {
      if($translation_count > $stop_codons_count) {
        StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::TRANSLATES |
                     StatMsg::PARTIAL);
      }
    } else {
      StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::TRANSLATES |
                   StatMsg::ENTIRE);
    }
  }
}


=head2 transfer_transcript

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub transfer_transcript {
    my $transcript = shift;
    my $mapper = shift;
    my $V_cs = shift;
    my $V_pfa = shift;
    my $E_slice = shift;

    $support->log_verbose("Transcript: " . $transcript->stable_id."\n", 3);

    my $V_exons = $transcript->get_all_Exons;
    my $E_cdna_pos = 0;
    my $cdna_exon_start = 1;

    my $E_transcript = InterimTranscript->new;
    $E_transcript->stable_id($transcript->stable_id);
    $E_transcript->version($transcript->version);
    $E_transcript->biotype($transcript->biotype);
    $E_transcript->confidence($transcript->confidence);
    $E_transcript->description($transcript->description);
    $E_transcript->created_date($transcript->created_date);
    $E_transcript->modified_date($transcript->modified_date);
    $E_transcript->cdna_coding_start($transcript->cdna_coding_start);
    $E_transcript->cdna_coding_end($transcript->cdna_coding_end);

    # protein features
    if (defined($transcript->translation)) {
        $E_transcript->add_ProteinFeatures(@{ $V_pfa->fetch_all_by_translation_id($transcript->translation->dbID) });
    }

    my @E_exons;

    EXON:
    foreach my $V_exon (@{ $V_exons }) {
        $support->log_verbose("Exon: " . $V_exon->stable_id . " chr=" . 
                $V_exon->slice->seq_region_name . " start=". 
                $V_exon->seq_region_start."\n", 4);

        my $E_exon = InterimExon->new;
        $E_exon->stable_id($V_exon->stable_id);
        $E_exon->version($V_exon->version);
        $E_exon->cdna_start($cdna_exon_start);
        $E_exon->start_phase($V_exon->phase);
        $E_exon->end_phase($V_exon->end_phase);

        # supporting evidence
        foreach my $sf (@{ $V_exon->get_all_supporting_features }) {
            # map coordinates
            my @coords = $mapper->map(
                    $sf->seq_region_name,
                    $sf->seq_region_start,
                    $sf->seq_region_end,
                    $sf->seq_region_strand,
                    $V_cs,
            );
            if (@coords == 1) {
                my $c = $coords[0];
                unless ($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
                    $sf->start($c->start);
                    $sf->end($c->end);
                    $sf->strand($c->strand);
                    $sf->slice($E_slice);
                    $E_exon->add_SupportingFeature($sf);
                }
            }
        }

        # map exon coordinates
        my @coords = $mapper->map(
                $V_exon->seq_region_name,
                $V_exon->seq_region_start,
                $V_exon->seq_region_end,
                $V_exon->seq_region_strand,
                $V_cs,
        );

        if (@coords == 1) {
            my $c = $coords[0];

            if ($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
            #
            # Case 1: Complete failure to map exon
            #

                $E_exon->fail(1);
                $E_transcript->add_Exon($E_exon);

            } else {
            #
            # Case 2: Exon mapped perfectly
            #

                $E_exon->start($c->start);
                $E_exon->end($c->end);
                $E_exon->cdna_start($cdna_exon_start);
                $E_exon->cdna_end($cdna_exon_start + $E_exon->length - 1);
                $E_exon->strand($c->strand);
                $E_exon->seq_region($c->id);

                $E_cdna_pos += $c->length;
            }
        } else {
        #
        # Case 3 : Exon mapped partially
        #
            
            $E_exon->fail(1);
            $E_transcript->add_Exon($E_exon);
        }

        $cdna_exon_start = $E_cdna_pos + 1;

        $E_transcript->add_Exon($E_exon);
    } # foreach exon

    return $E_transcript;
}


=head2 create_transcripts

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub create_transcripts {
    my $itranscript   = shift;
    my $E_sa = shift;

    # check the exons and split transcripts where exons are bad
    my $itranscripts = Transcript::check_iexons($support, $itranscript);

    my @finished_transcripts;
    my %protein_features;
    foreach my $itrans (@{ $itranscripts }) {
        # if there are any exons left in this transcript add it to the list
        if (@{ $itrans->get_all_Exons }) {
            my ($tr, $pf) = Transcript::make_Transcript($support, $itrans,
                $E_sa);
            push @finished_transcripts, $tr;
            $protein_features{$tr->stable_id} = $pf;
        } else {
            $support->log("Transcript ". $itrans->stable_id . " has no exons left.\n", 3);
        }
    }

    return \@finished_transcripts, \%protein_features;
}


