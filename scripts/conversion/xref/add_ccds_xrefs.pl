#!/usr/local/bin/perl

=head1 NAME

add_ccds_xrefs.pl - adds xrefs for CCDS identifiers

=head1 SYNOPSIS

add_ccds_xrefs.pl [options]

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
    --ccdsfile FILE                     read input from FILE

=head1 DESCRIPTION

This script adds CCDS identifiers to the database. The input file is a
whitespace-separated list of transcript stable IDs and CCDS identifiers. For
more information on the CCDS database see
http://www.ensembl.org/Homo_sapiens/ccds.html.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHORS

Steve Trevanion <st3@sanger.ac.uk>
Patrick Meidl <pm2@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);
use Storable;
use Data::Dumper;

BEGIN {
  $SERVERROOT = "$Bin/../../../..";
  unshift(@INC, "$SERVERROOT/ensembl/modules");
  unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
  'prune=s',
  'evegahost=s',
  'evegaport=s',
  'evegauser=s',
  'evegapass=s',
  'evegassembly=s',
  'ensemblhost=s',
  'ensemblport=s',
  'ensembluser=s',
  'ensemblpass=s',
  'ccdshost=s',
  'ccdsport=s',
  'ccdsuser=s',
  'ccdspass=s',
);
$support->allowed_params(
  $support->get_common_params,
  'prune',
  'evegahost',
  'evegaport',
  'evegauser',
  'evegapass',
  'evegaassembly',
  'evegadbname',
  'ensemblhost',
  'ensemblport',
  'ensembluser',
  'ensemblpass',
  'ccdshost',
  'ccdsport',
  'ccdsuser',
  'ccdspass',
  'ccdsdbname',
);

if ($support->param('help') or $support->error) {
  warn $support->error if $support->error;
  pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# get adaptors for adding to Vega
my $dba = $support->get_database('ensembl');
my $ta = $dba->get_TranscriptAdaptor;
my $ea = $dba->get_DBEntryAdaptor;

if ($support->param('prune') and $support->user_proceed('Would you really like to delete all previous CCDS xrefs?')) {
  $support->log("Deleting  CCDS  xrefs...\n");
  my $num = $dba->dbc->do(qq(
           DELETE x
           FROM xref x, external_db ed
           WHERE x.external_db_id = ed.external_db_id
           AND ed.db_name = 'CCDS'));
  $support->log("Done deleting $num entries.\n");

  # object_xrefs
  $support->log("Deleting orphan object_xrefs...\n");
  $num = $dba->dbc->do(qq(
           DELETE ox
           FROM object_xref ox
           LEFT JOIN xref x ON ox.xref_id = x.xref_id
           WHERE x.xref_id IS NULL
        ));
}

#parse ccds info
my ($vega_ids,$ccds_ids);
my $vega_ids_file = $support->param('logpath').'/'.$support->param('ccdsdbname')."-parsed_vega_ids_for_ccds.file";
if (-e $vega_ids_file) {
  if ($support->user_proceed("Read CCDS records from a previously saved file ($vega_ids_file) ?\n")) {
    $vega_ids = retrieve($vega_ids_file);
  }
}
else {
  my $e_dba = $support->get_database('ensembl','ensembl');
  my $ev_dba = $support->get_database('ensembl','evega');
  my $ev_sa = $ev_dba->get_SliceAdaptor;
  my $ccds_dba = $support->get_database('ensembl','ccds');
  my $ccds_sa = $ccds_dba->get_SliceAdaptor;
  $ev_dba->dnadb($e_dba);
  $ccds_dba->dnadb($e_dba);
  my $path = $support->param('evegaassembly');
  $support->log_stamped("Retrieving info from CCDS database\n");
  foreach my $chr ( @{$ccds_sa->fetch_all('chromosome')}){
    my $chrname = $chr->name;
    $support->log_stamped("On chromosome $chrname\n",1);
    foreach my $ccds_gene ( @{$chr->get_all_Genes()} ){
      $ccds_gene = $ccds_gene->transform('chromosome', $path);
      foreach my $ccds_trans (@{$ccds_gene->get_all_Transcripts()}){
        my %xref_hash;
        my $ccds_id;
        foreach my $entry (@{$ccds_trans->get_all_DBEntries('CCDS')}) {
          $xref_hash{$entry->display_id()} = 1;
        }
        if (scalar keys %xref_hash != 1){
          my $tsi = $ccds_trans->stable_id;
          $support->log_error("Something odd going on for $tsi, has ". scalar (keys %xref_hash) ." xrefs\n",1);
          foreach my $entry (keys %xref_hash){
            $support->log("xref $entry \n",2);
          }
        }
        else {
          foreach my $entry (keys %xref_hash){
            $ccds_id = $entry;
          }
        }
        my $chr_name = $ccds_trans->slice->seq_region_name;
        my $start = $ccds_trans->start();
        my $end   = $ccds_trans->end();
        my @ccds_exons = @{$ccds_trans->get_all_translateable_Exons()};
        my $slice = $ev_sa->fetch_by_region('chromosome',$chr_name, $start, $end, '1', $path );
        foreach my $gene (@{$slice->get_all_Genes()}){
          next if ($gene->biotype ne "protein_coding");
          $gene = $gene->transform('chromosome', $path);
          foreach my $trans (@{$gene->get_all_Transcripts}){
            my @exons = @{$trans->get_all_translateable_Exons()};
            my $match = 0;
            if (scalar @exons == scalar @ccds_exons){
              for (my $i = 0; $i < @exons; $i++){
                if ($ccds_exons[$i]->start == $exons[$i]->start
                      && $ccds_exons[$i]->end == $exons[$i]->end){
                  $match++;
                }
                #else{
                #  print "no match ".$ccds_exons[$i]->start." != ".$exons[$i]->start." or ".
                #	$ccds_exons[$i]->end." != ".$exons[$i]->end."\n";
                #}
              }
              if ($match == scalar @exons){
                $vega_ids->{$trans->stable_id} = $ccds_id;
              }
            }
          }
        }
      }
    }
  }
  store($vega_ids,$vega_ids_file);
}

foreach my $vid (keys %$vega_ids) {
  $ccds_ids->{$vega_ids->{$vid}}++;
}
$support->log_stamped("Done parsing, retrieved ".scalar(keys %$vega_ids)." transcripts and ".scalar(keys %$ccds_ids)." CCDS entries\n");

# loop over Vega transcripts, updating db with CCDS identifier
$support->log_stamped("Adding xrefs to vega db...\n");
my ($no_trans, $num_success, $no_match) = (0, 0, 0);
my (%sources, %transcript_ids);
my ($successful, $non_translating, $missing_transcript);
while (my ($tsi, $ccds_idnt) = each %$vega_ids) {
  if (my $transcript = $ta->fetch_by_stable_id($tsi)) {
    if (my $translation = $transcript->translation) {
      $sources{$transcript->analysis->logic_name}++;
      $transcript_ids{$tsi}++;
      $num_success++;
      my $internal_id = $translation->dbID;
      my ($prim_acc) = $ccds_idnt =~ /(\w*)/;
      $successful .= sprintf "    %-30s%-20s%-20s\n", $tsi, $internal_id, $ccds_idnt;
      my $dbentry = Bio::EnsEMBL::DBEntry->new(
	-primary_id => $prim_acc,
	-display_id => $ccds_idnt,
	-version    => 1,
	-dbname     => 'CCDS',
      );
      unless ($support->param('dry_run')) {
	$ea->store($dbentry, $internal_id, 'Translation');
      }
    } else {
      $no_trans++;
      $non_translating .= "    $tsi\n";
    }
  } else {
    $no_match++;
    $missing_transcript .= "    $tsi\n";
  }
}
$support->log("Done. ".$support->date_and_mem."\n\n");

# print log results
$support->log("\nProcessed ".scalar(keys %$vega_ids)." identifiers.\n");
$support->log("OK: $num_success\n");
$support->log("WARNINGS:\n");
$support->log("Identifiers with no matching transcript in Vega: $no_match.\n", 1);
$support->log("Transcripts in this set that don't translate: $no_trans.\n", 1);

if ($successful) {
  $support->log("\nTranscripts which had a CCDS identifier added:\n");
  $support->log(sprintf "    %-30s%-20s%-20s\n", qw(STABLE_ID DBID CCDS_ID));
  $support->log("    " . "-"x70 . "\n");
  $support->log($successful);
}
if ($missing_transcript) {
  $support->log_warning("\nTranscripts with no matching CCDS identifier:\n");
  $support->log($missing_transcript);
}
if ($non_translating) {
  $support->log_warning("\nTranscripts in this set that don't translate:\n");
  $support->log($non_translating);
}
$support->log("\nTranscripts in with xrefs added had the following logic_names:\n");
while (my ($ln, $c) = each %sources) {
  $support->log("$c with logic_name of $ln\n",1);
}
my @dup_ids = grep {$transcript_ids{$_} > 1} keys %transcript_ids;
if (@dup_ids) {
  $support->log("There are scalar(@dup_ids) with duplicates:\n" . join "\n" , @dup_ids . "\n");
}


# finish log
$support->finish_log;

