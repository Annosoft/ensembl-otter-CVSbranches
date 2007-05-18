#!/usr/local/bin/perl

=head1 NAME

examine_GD.pl - identify redundant GD genes for deleting

=head1 SYNOPSIS

examine_GD.pl [options]

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

    --delete=FILE                  file for stable IDs of GD genes to delete
    --chromosomes=LIST                  list of chromosomes to read (not working)

=head1 DESCRIPTION

examine GD genes to identify redundant ones (ie have an overlapping Havana gene) and generate
a file of stable IDs to be used for deleting. Verbosely report numbers of GD genes with a
gomi remark (for jla1).

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
use FindBin qw($Bin);
use vars qw($SERVERROOT);

#use Data::Dumper;

BEGIN {
    $SERVERROOT = "$Bin/../../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Analysis;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'chromosomes|chr=s@',
	'delete=s',
);
$support->allowed_params(
    $support->get_common_params,
    'chromosomes',
	'delete',
);
if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');
$support->param('delete',($support->param('logpath').'/'.$support->param('delete')));

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# connect to database and get adaptors
my $dba = $support->get_database('otter');
my $dbh = $dba->dbc->db_handle;
my $ga  = $dba->get_GeneAdaptor;
my $sa  = $dba->get_SliceAdaptor;
my $ta  = $dba->get_TranscriptAdaptor;

#filehandle for output
my $outfile = $support->filehandle('>', $support->param('delete'));

#hashes for more detailed logging if ever needed
my (%to_delete,%gomi_to_log_overlap,%gomi_to_log_no_overlap);
#counters
my ($tot_c,$noverlap_c,$overlap_c);

CHR:
foreach my $chr (@{$sa->fetch_all('chromosome')}) {
	next CHR if ($chr->seq_region_name eq 'chrY-03');
	$support->log_stamped("Looping over chromosome ".$chr->seq_region_name."\n");
 GENE:
    foreach my $gene (@{ $ga->fetch_by_Slice($chr) }) {
		my $gsi = $gene->stable_id;
		my $name = $gene->gene_info->name->name;
		next GENE if ($name !~ /^GD:/);
		$tot_c++;
		$support->log_verbose("Studying gene $gsi ($name)\n");
	
		#get any overlapping Havana genes
		my $slice = $gene->feature_Slice;
		if (my @genes = @{$ga->fetch_all_by_Slice($slice,'otter')}) {
			$overlap_c++;
			#note stable id of file to be deleted
			$to_delete{$gsi} = $name;
			#are there any gomi remarks ?
			foreach my $trans (@{$gene->get_all_Transcripts()}) {
				my $tsi = $trans->stable_id;
				if (&any_gomi_remarks($trans)) {
					push @{$gomi_to_log_overlap{$gsi}}, $tsi;
				}
			}
		}
		else {
			$noverlap_c++;
			#are there any gomi remarks ?
			foreach my $trans (@{$gene->get_all_Transcripts()}) {
				my $tsi = $trans->stable_id;
				if (&any_gomi_remarks($trans)) {
					push @{$gomi_to_log_no_overlap{$gsi}}, $tsi;
				}
			}
		}
	}			
}

#logging
my $log_c_o = keys %gomi_to_log_overlap;
my $log_c_no = keys %gomi_to_log_no_overlap;
$support->log("There are $tot_c GD genes in total:\n");	
$support->log("$noverlap_c of these do not overlap with Havana\n",1);
$support->log("$overlap_c of these overlap with Havana and wil be deleted\n",1);

$support->log_verbose("$log_c_no GD genes do not overlap Havana loci but do have a gomi remark\n");
$support->log_verbose("$log_c_o GD genes overlap Havana loci (and will be pruned from Vega) and have a gomi remark\n");

#create file to delete
print $outfile join("\n", keys %to_delete), "\n";
close $outfile;

$support->finish_log;

sub any_gomi_remarks {
	my ($trans) = @_;
	my $info = $trans->transcript_info;
	foreach my $remark ($info->remark) {
		my $value = $remark->remark;
		return 1 if ($value =~ /gomi/);	
	}
}
