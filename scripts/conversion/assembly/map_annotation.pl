#!/usr/local/bin/perl

=head1 NAME

map_annotation.pl - map features from one assembly onto another

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
    --chromosomes, --chr=LIST           only process LIST chromosomes
    --prune=0|1                         delete results from previous runs of
                                        this script first
    --logic_names=LIST                  restrict transfer to gene logic_names

=head1 DESCRIPTION

This script is part of a series of scripts to transfer annotation from a
Vega to an Ensembl assembly. See "Related scripts" below for an overview of the
whole process.

Given a database with a mapping between two different assemblies of the same
genome, this script transfers features from one assembly to the other.

Features transfer include:

    - genes/transcripts/exons/translations
    - xrefs
    - supporting features and associated dna/protein_align_features
    - protein features

Currently, only complete transfers are considered. This is the easiest way to
ensure that the resulting gene structures are identical to the original ones.
For future release, there are plans to store incomplete matches by using the
Ensembl API's SeqEdit facilities.

Genes transferred can be restricted on logic_name using the --logic_names
option. Used for mouse (-logic_names otter,otter_external).

Look in the logs for 'Set coordinates' and check exon coordinates of any examples
- untested code.

=head1 RELATED SCRIPTS

The whole Ensembl-vega database production process is done by these scripts:

    ensembl-otter/scripts/conversion/assembly/make_ensembl_vega_db.pl
    ensembl-otter/scripts/conversion/assembly/map_annotation.pl
    ensembl-otter/scripts/conversion/assembly/finish_ensembl_vega_db.pl

See documention in the respective script for more information.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>
Based on code originally wrote by Graham McVicker and Patrick Meidl

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
use Deletion;
use Transcript;
use Gene;

$| = 1;

our $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

#$SIG{INT}= 'do_stats_logging';      # signal handler for Ctrl-C, i.e. will call sub do_stats_logging


# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
  'evegahost=s',
  'evegaport=s',
  'evegauser=s',
  'evegapass=s',
  'evegadbname=s',
  'chromosomes|chr=s@',
  'logic_names=s@',
  'prune=s',
);
$support->allowed_params(
  $support->get_common_params,
  'evegahost',
  'evegaport',
  'evegauser',
  'evegapass',
  'evegadbname',
  'chromosomes',
  'logic_names',
  'prune',
);

if ($support->param('help') or $support->error) {
  warn $support->error if $support->error;
  pod2usage(1);
}

$support->comma_to_list('chromosomes');
$support->comma_to_list('logic_names');
my @logic_names = ( $support->param('logic_names'));

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

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
my $E_ta = $E_dba->get_TranscriptAdaptor;
my $E_pfa = $E_dba->get_ProteinFeatureAdaptor;
my $cs_adaptor = $E_dba->get_CoordSystemAdaptor;
my $asmap_adaptor = $E_dba->get_AssemblyMapperAdaptor;

my $E_cs = $cs_adaptor->fetch_by_name('chromosome',
    $support->param('ensemblassembly'));
my $V_cs = $cs_adaptor->fetch_by_name('chromosome',
    $support->param('assembly'));

# get assembly mapper
my $mapper = $asmap_adaptor->fetch_by_CoordSystems($E_cs, $V_cs);
$mapper->max_pair_count( 6_000_000 );
$mapper->register_all;

# if desired, delete entries from previous runs of this script
if ($support->param('prune') && $support->user_proceed("Do you want to delete all entries from previous runs of this script?")) {
  $support->log("Deleting db entries from previous runs of this script...\n");
  $E_dbh->do(qq(DELETE FROM analysis));
  $E_dbh->do(qq(DELETE FROM dna_align_feature));
  $E_dbh->do(qq(DELETE FROM exon));
  $E_dbh->do(qq(DELETE FROM exon_stable_id));
  $E_dbh->do(qq(DELETE FROM exon_transcript));
  $E_dbh->do(qq(DELETE FROM gene));
  $E_dbh->do(qq(DELETE FROM gene_stable_id));
  $E_dbh->do(qq(DELETE FROM gene_attrib));
  $E_dbh->do(qq(DELETE FROM object_xref));
  $E_dbh->do(qq(DELETE FROM protein_align_feature));
  $E_dbh->do(qq(DELETE FROM protein_feature));
  $E_dbh->do(qq(DELETE FROM transcript_supporting_feature));
  $E_dbh->do(qq(DELETE FROM supporting_feature));
  $E_dbh->do(qq(DELETE FROM transcript));
  $E_dbh->do(qq(DELETE FROM transcript_stable_id));
  $E_dbh->do(qq(DELETE FROM transcript_attrib));
  $E_dbh->do(qq(DELETE FROM translation));
  $E_dbh->do(qq(DELETE FROM translation_stable_id));
  $E_dbh->do(qq(DELETE FROM translation_attrib));
  $E_dbh->do(qq(DELETE x
                  FROM xref x, external_db ed
                  WHERE x.external_db_id = ed.external_db_id
                  AND ed.db_name NOT IN ('Interpro')
     ));
  $support->log("Done.\n");
}

my (%stat_hash,%trans_numbers);

# loop over chromosomes
$support->log("Looping over chromosomes...\n");
my $V_chrlength = $support->get_chrlength($E_dba, $support->param('assembly'),'chromosome',1);
my $E_chrlength = $support->get_chrlength($E_dba, $support->param('ensemblassembly'),'chromosome',1);
my $ensembl_chr_map = $support->get_ensembl_chr_mapping($V_dba, $support->param('assembly'));

foreach my $V_chr ($support->sort_chromosomes($V_chrlength)) {

  $support->log_stamped("Chromosome $V_chr...\n", 1);

  # skip non-ensembl chromosomes
  my $E_chr = $ensembl_chr_map->{$V_chr};
  unless ($E_chrlength->{$E_chr}) {
    $support->log_warning("Ensembl chromosome equivalent to $V_chr not found. Skipping.\n", 1);
    next;
  }

  # fetch chromosome slices
  my $V_slice = $V_sa->fetch_by_region('chromosome', $V_chr, undef, undef,
				       undef, $support->param('assembly'));
  my $E_slice = $E_sa->fetch_by_region('chromosome', $E_chr, undef, undef,
				       undef, $support->param('ensemblassembly'));

  $support->log("Looping over genes...\n", 1);
  my $genes = $V_ga->fetch_all_by_Slice($V_slice);
 GENE:
  foreach my $gene (@{ $genes }) {
    my $gsi = $gene->stable_id;
    my $ln = $gene->analysis->logic_name;
    my $name = $gene->display_xref->display_id;
    if ($support->param('logic_names')) {
      unless (grep {$ln eq $_} @logic_names) {
	$support->log_verbose("Skipping gene $gsi/$name (logic_name $ln)\n",2);
	next GENE;
      }
    }
    $support->log("Gene $gsi/$name (logic_name $ln)\n", 2);

    my $transcripts = $gene->get_all_Transcripts;
    my (@finished, %all_protein_features);
    my $c = 0;
  TRANS:
    foreach my $transcript (@{ $transcripts }) {
      if ($transcript->biotype eq 'artifact') {
	$support->log("Transcript: ".$transcript->stable_id." skipping because it's an 'artifact'\n", 3);
	next TRANS;
      }
      $c++;

      my $interim_transcript = transfer_transcript($transcript, $mapper,
						   $V_cs, $V_pfa, $E_slice);
      my ($finished_transcripts, $protein_features) =
	create_transcripts($interim_transcript, $E_sa, $gsi);

      # set the translation stable identifier on the finished transcripts
      foreach my $tr (@{ $finished_transcripts }) {
	if ($tr->translation && $transcript->translation) {
	  $tr->translation->stable_id($transcript->translation->stable_id);
	  $tr->translation->version($transcript->translation->version);
	  $tr->translation->created_date($transcript->translation->created_date);
	  $tr->translation->modified_date($transcript->translation->modified_date);
	}
      }
			
      push @finished, @$finished_transcripts;
      map { $all_protein_features{$_} = $protein_features->{$_} }
	keys %{ $protein_features || {} };
    }

    # if there are no finished transcripts, count this gene as being NOT transfered
    my $num_finished_t= @finished;
    if(! $num_finished_t){
      push @{$stat_hash{$V_chr}->{'failed'}}, [$gene->stable_id,$gene->seq_region_start,$gene->seq_region_end];
      next GENE;
    }

    #make a note of the number of transcripts per gene
    $trans_numbers{$gsi}->{'vega'} = scalar(@{$transcripts});
    $trans_numbers{$gsi}->{'evega'} = $num_finished_t;

    #count gene and transcript if it's been transferred
    $stat_hash{$V_chr}->{'genes'}++;
    $stat_hash{$V_chr}->{'transcripts'} += $c;
	
    unless ($support->param('dry_run')) {
      Gene::store_gene($support, $E_slice, $E_ga, $E_ta, $E_pfa, $gene,
		       \@finished, \%all_protein_features);
    }
  }
  $support->log("Done with chromosome $V_chr.\n", 1);
}

#see if any transcripts / gene are different
foreach my $gsi (keys %trans_numbers) {
  if ($trans_numbers{$gsi}->{'vega'} != $trans_numbers{$gsi}->{'evega'}) {
    my $v_num = $trans_numbers{$gsi}->{'vega'};
    my $e_num = $trans_numbers{$gsi}->{'evega'};
    $support->log("There are different numbers of transcripts for gene $gsi in vega ($v_num) and ensembl_vega ($e_num)\n");
  }
}

# write out to statslog file
do_stats_logging();

# finish logfile
$support->finish_log;

### END main


=head2 transfer_transcript

  Arg[1]      : Bio::EnsEMBL::Transcript $transcript - Vega source transcript
  Arg[2]      : Bio::EnsEMBL::ChainedAssemblyMapper $mapper - assembly mapper
  Arg[3]      : Bio::EnsEMBL::CoordSystem $V_cs - Vega coordinate system
  Arg[4]      : Bio::EnsEMBL::ProteinFeatureAdaptor $V_pfa - Vega protein
                feature adaptor
  Arg[5]      : Bio::EnsEMBL::Slice $slice - Ensembl slice
  Description : This subroutine takes a Vega transcript and transfers it (and
                all associated features) to the Ensembl assembly.
  Return type : InterimTranscript - the interim transcript object representing
                the transfered transcript
  Exceptions  : none
  Caller      : internal

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
  $E_transcript->status($transcript->status);
  $E_transcript->description($transcript->description);
  $E_transcript->created_date($transcript->created_date);
  $E_transcript->modified_date($transcript->modified_date);
  $E_transcript->cdna_coding_start($transcript->cdna_coding_start);
  $E_transcript->cdna_coding_end($transcript->cdna_coding_end);
  $E_transcript->transcript_attribs($transcript->get_all_Attributes);
  $E_transcript->analysis($transcript->analysis);

  # transcript supporting evidence
  foreach my $sf (@{ $transcript->get_all_supporting_features }) {
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
	$E_transcript->add_TranscriptSupportingFeature($sf);
      }
    }
  }

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
    $E_exon->created_date($V_exon->created_date);
    $E_exon->modified_date($V_exon->modified_date);
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

      }
      else {
	#
	# Case 2: Exon mapped perfectly
	#
	$E_exon->start($c->start);
	$E_exon->end($c->end);
	$E_exon->strand($c->strand);
	$E_exon->seq_region($c->id);

	$E_exon->cdna_start($cdna_exon_start);
	$E_exon->cdna_end($cdna_exon_start + $E_exon->length - 1);
	$E_cdna_pos += $c->length;
      }
    }
    else {
      my ($start,$end,$strand,$seq_region);
      my $gap = 0;	
      foreach my $c (@coords) {
      #
      # Case 3 : Exon mapped partially
      #
	if ($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
	  $E_exon->fail(1);
	  $support->log_warning("Reason: Exon ".$E_exon->stable_id." mapping has a gap\n",4);
	  $gap = 1;
	  last;
	}
      }
      #
      # Case 4 : Multiple mappings for exon, but no gaps
      #
      unless ($gap) {
	my ($last_end);
	foreach my $c (sort {$a->start <=> $b->start} @coords) {
	  if ($last_end) {
	    if ($c->start != $last_end) {
	      $E_exon->fail(1);
	      $support->log_warning("Reason: Exon mapping has a mismatch in coords\n",4);
	      last;
	    }
	  }
	  $start = ! $start ? $c->start
	           : $start > $c->start ? $c->start
	           : $start;
	  $end = ! $end ? $c->end
	         : $end < $c->end ? $c->end
	         : $end;
	  $strand = $c->strand;
	  $seq_region = $c->id;
	  $last_end = $c->end;
	}
	
	$E_exon->start($start);
	$E_exon->end($end);
	$E_exon->strand($strand);
	$E_exon->seq_region($seq_region);

	$E_exon->cdna_start($cdna_exon_start);	
	my $length = $end-$start+1;
	$E_exon->cdna_end($cdna_exon_start + $length);
	$E_cdna_pos += $length+1;

	unless ($E_exon->fail) {
	  $support->log_warning("Set coordinates but this is untested code, please check\n");
	}

      }
    }
    $cdna_exon_start = $E_cdna_pos + 1;

    $E_transcript->add_Exon($E_exon);
  }
  return $E_transcript;
}



=head2 create_transcripts

  Arg[1]      : InterimTranscript $itranscript - an interim transcript object
  Arg[2]      : Bio::EnsEMBL::SliceAdaptor $E_sa - Ensembl slice adaptor
  Arg[3]      : Gene Stable ID (used for checking against curation)
  Description : Creates the actual transcripts from interim transcripts
  Return type : List of a listref of Bio::EnsEMBL::Transcripts and a hashref of
                protein features (keys: transcript_stable_id, values:
                Bio::Ensmembl::ProteinFeature)
  Exceptions  : none
  Caller      : internal

=cut

sub create_transcripts {
  my $itranscript   = shift;
  my $E_sa = shift;
  my $gsi = shift;

  # check the exons and split transcripts where exons are bad
  my $itranscripts = Transcript::check_iexons($support, $itranscript, $gsi);

  my @finished_transcripts;
  my %protein_features;
  foreach my $itrans (@{ $itranscripts }) {
    # if there are any exons left in this transcript add it to the list
    if (@{ $itrans->get_all_Exons }) {
      my ($tr, $pf) = Transcript::make_Transcript($support, $itrans, $E_sa);
      push @finished_transcripts, $tr;
      $protein_features{$tr->stable_id} = $pf;
    } else {
      $support->log("Transcript ". $itrans->stable_id . " has no exons left.\n", 3);
    }
  }
  return \@finished_transcripts, \%protein_features;
}

sub do_stats_logging{

  #writes the number of genes and transcripts processed to the log file
  #note: this can be called as an interrupt handler for ctrl-c,
  #so can also give current stats if script terminated

  my %failed;
  my $format = "%-20s%-10s%-10s\n";
  $support->log(sprintf($format,'Chromosome','Genes','Transcripts'));
  my $sep = '-'x41;
  $support->log("$sep\n");
  foreach my $chrom(sort keys %stat_hash){
    my $num_genes= $stat_hash{$chrom}->{'genes'};
    my $num_transcripts= $stat_hash{$chrom}->{'transcripts'};
    if(defined($stat_hash{$chrom}->{'failed'})){
      $failed{$chrom} = $stat_hash{$chrom}->{'failed'};
    }
    $support->log(sprintf($format,$chrom,$num_genes,$num_transcripts));
  }
  $support->log("\n");
  foreach my $failed_chr (keys %failed){
    my $no = scalar @{$failed{$failed_chr}};
    $support->log("$no genes not transferred on chromosome $failed_chr:\n");
    foreach my $g (@{$failed{$failed_chr}}) {
      $support->log("  ".$g->[0].": ".$g->[1]."-".$g->[2]."\n");
    }
  }
}
