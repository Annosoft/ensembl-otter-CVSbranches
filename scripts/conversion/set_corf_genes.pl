#!/usr/local/bin/perl

#!/usr/local/bin/perl

=head1 NAME

set_corf_genes.pl - set analysis_id for CORF genes and transcripts

=head1 SYNOPSIS

set_corf_genes.pl [options]

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

    --logic_name=NAME                   logicname for CORF genes (defaults to otter_corf)
    --chromosomes=LIST                  List of chromosomes to read (not working)

=head1 DESCRIPTION

This script identifies all CORF genes in a Vega database and sets the analysis_id for each.
It also sets the analysis_id of the transcripts.

Logic is:
(i) retrieve each gene with a GD: prefix on it's name (gene_name.name)
(ii) if that gene has an 'Annotation_remark- corf' remark on at least one of it's transcripts
then set the analysis_id to otter_corf
(iii) if that gene has a transcript with a 'corf' remark on it then change that to be
'Annotation_remark- corf', and then set the analysis_id as above

Reports on other GD: loci that have other remarks containing 'corf'.
Verbosely reports on non GD: loci that have remarks containing 'corf'.

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
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);

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

our $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'logic_name=s',
    'chromosomes|chr=s@',
);
$support->allowed_params(
    $support->get_common_params,
    'logic_name',
    'chromosomes',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');

$support->param('logic_name','otter_corf') unless $support->param('logic_name');

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

my $sth1 = $dbh->prepare("insert into transcript_remark values ('',?,?)");
my $sth2 = $dbh->prepare("delete from transcript_remark where transcript_info_id = ? and remark = ?");
my (%non_GD,%to_change,%wrong_syntax_GD,%wrong_syntax_nonGD);
CHR:
foreach my $chr (@{$sa->fetch_all('chromosome')}) {
#	next CHR if ($chr->seq_region_name eq 'chrY-03');
	$support->log_stamped("Looping over chromosome ".$chr->seq_region_name."\n");
 GENE:
    foreach my $gene (@{ $ga->fetch_by_Slice($chr) }) {
		my $gsi = $gene->stable_id;
		my $name = $gene->gene_info->name->name;
		$support->log_verbose("Studying gene $gsi ($name)\n");
		my $comment = 0;
		foreach my $trans (@{$gene->get_all_Transcripts()}) {
			my $tsi = $trans->stable_id;
			my $info = $trans->transcript_info;
			my $t_info_id = $info->dbID;
			foreach my $remark ($info->remark) {
				my $value = $remark->remark;
				if ($value =~ /corf/i) {
					if ( $value =~ /^Annotation_remark- corf/i) {
						#capture nonGD genes with the correct remark
						if ($name !~ /^GD:/) {
							$non_GD{$gsi}->{'name'} = $name;
							push @{$non_GD{$gsi}->{'transcripts'}},$tsi;
						}
						#capture genes to patch analysis_id
						else {							
							$to_change{$gsi}->{'name'} = $name;
							push @{$to_change{$gsi}->{'transcripts'}}, [$tsi,$value];
						}
					}
					#capture genes with a remark that should be a corf annotation remark
					elsif ( ($value =~ /^corf/i) && ($name =~ /^GD:/)) { 
						unless (! $support->param('dry_run')) {
							$sth1->execute('Annotation_remark- corf',$t_info_id) unless (! $support->param('dry_run'));
							$sth2->execute($t_info_id,'corf');
						}
						$support->log_warning("Setting gene $gsi to corf despite not having the properly formatted Annotation remark\n");
						$to_change{$gsi}->{'name'} = $name;
						push @{$to_change{$gsi}->{'transcripts'}}, [$tsi,$value];
					}
					#capture incorrect syntax
					else {
						if ($name =~ /^GD:/) {
							$wrong_syntax_GD{$gsi}->{'name'} = $name;
							push @{$wrong_syntax_GD{$gsi}->{'transcripts'}}, [$tsi,$value];
						}
						else {
							$wrong_syntax_nonGD{$gsi}->{'name'} = $name;
							push @{$wrong_syntax_nonGD{$gsi}->{'transcripts'}}, [$tsi,$value];
						}
					}
				}
			}									
		}
	}
}

#report on GD genes with the incorrect syntax
$support->log("\nThe following are GD genes with the incorrect syntax:\n"); 
foreach my $gsi (keys %wrong_syntax_GD) {
	$support->log("\n$gsi (".$wrong_syntax_GD{$gsi}->{'name'}."):\n",1);
	foreach my $t (@{$wrong_syntax_GD{$gsi}->{'transcripts'}}) {
		$support->log($t->[0]." (".$t->[1].")\n",2);
	}
}

#report on non GD genes with the incorrect syntax
$support->log_verbose("\nThe following are non GD genes with the incorrect CORF syntax:\n"); 
foreach my $gsi (keys %wrong_syntax_nonGD) {
	$support->log_verbose("\n$gsi (".$wrong_syntax_nonGD{$gsi}->{'name'}."):\n",1);
	foreach my $t (@{$wrong_syntax_nonGD{$gsi}->{'transcripts'}}) {
		$support->log_verbose($t->[0]." (".$t->[1].")\n",2);
	}
}

#report on non GD genes with the correct syntax
$support->log_verbose("\nThe following are non GD genes with the correct CORF syntax:\n");
my $c = 0;
foreach my $gsi (keys %non_GD) {
	$c++;
	$support->log_verbose("\n$gsi (".$non_GD{$gsi}->{'name'}."):\n",1);
	foreach my $t (@{$non_GD{$gsi}->{'transcripts'}}) {
		$support->log_verbose($t."\n",2);
	}
}
$support->log("\n$c non GD genes with the correct CORF syntax in total.\n\n");

#report on GD genes with the correct syntax, ie the ones to be updated
$support->log("\nThe following are GD genes with the correct CORF syntax and will be updated:\n");
$c = 0;
my @gene_stable_ids;
foreach my $gsi (keys %to_change) {
	$c++;
	push @gene_stable_ids, $gsi;
	$support->log("\n$gsi (".$to_change{$gsi}->{'name'}."):\n",1);
	foreach my $t (@{$to_change{$gsi}->{'transcripts'}}) {
		$support->log_verbose($t->[0]." (".$t->[1].")\n",2);
	}
}
$support->log("\n$c GD genes with the correct CORF syntax in total.\n\n");

#do the updates
# add analysis
if (! $support->param('dry_run')) {
    $support->log("Adding new analysis...\n");
    my $analysis = new Bio::EnsEMBL::Analysis (
        -program     => "set_corf_genes.pl",
        -logic_name  => $support->param('logic_name'),
    );
    my $analysis_id = $dba->get_AnalysisAdaptor->store($analysis);
    $support->log_error("Couldn't store analysis ".$support->param('analysis').".\n") unless $analysis_id;

    # change analysis for genes in list
    $support->log("Updating analysis of genes in list...\n");
    my $gsi_string = join("', '", @gene_stable_ids);
    my $num = $dbh->do(qq(
        UPDATE gene g, gene_stable_id gsi
        SET analysis_id = $analysis_id
        WHERE g.gene_id = gsi.gene_id
        AND gsi.stable_id in ('$gsi_string')
    ));
    $support->log("Done updating $num genes.\n\n");

	#change analysis_ids of transcripts
	if ($support->user_proceed("\nSet all transcript analysis_ids to equal those of their genes ?")) {	
		$support->log("Updating analysis of corresponding transcripts...\n");
		$dbh->do(qq(
            UPDATE transcript t, gene g
            SET t.analysis_id = g.analysis_id
            WHERE g.gene_id = t.gene_id
        ));
		$support->log("Done updating transcripts.\n\n");
	}
	else {
		$support->log("Transcripts analysis_ids not updated.\n\n");
	}
}

$support->finish_log;
