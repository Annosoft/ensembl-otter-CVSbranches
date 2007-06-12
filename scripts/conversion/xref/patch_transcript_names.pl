#!/usr/local/bin/perl

=head1 NAME

patch_transcript_names.pl - update transcript names with external names from genes

=head1 SYNOPSIS

patch_transcript_names.pl [options]

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

This script updates display xrefs for transcripts and translations to reflect changes
 in gene names since the xrefs were originally generated

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>

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
my $dba = $support->get_database('ensembl');
my $sa = $dba->get_SliceAdaptor;
my $ga = $dba->get_GeneAdaptor;

# statement handle for xref update
my $sth = $dba->dbc->prepare("update xref set display_label=? where xref_id=?");

# get genes
my $c;
foreach my $chr ($support->sort_chromosomes) {
	$support->log_stamped("Looping over chromosome $chr\n");
	my $slice = $sa->fetch_by_region('toplevel', $chr);
	foreach my $gene (@{$slice->get_all_Genes()}) {
		my $g_name = $gene->display_xref->display_id;
		my $gsi = $gene->stable_id;
		my $ln = $gene->analysis->logic_name;
		TRANS: foreach my $trans (@{$gene->get_all_Transcripts()}) {
			my $trans_dbentry = $trans->display_xref;
			my $stable_id =  $trans->stable_id;
			my ($t_name,$t_version) = $trans_dbentry->display_id =~ /(.*)-(\d+)$/;
			if (! $t_version) {
				if ($ln eq 'otter') {
					$support->log_warning("No correctly formatted version found for otter transcript $stable_id, please investigate. Not setting\n") unless ($t_version);
				}
				else {
					$support->log_verbose("WARNING: No correctly formatted version found for transcript $stable_id, not setting\n") unless ($t_version);
				}
				next TRANS;
			}
			if ($t_name ne $g_name) {
				$c++;
				my $new_name = $g_name.'-'.$t_version;
				$support->log("$c. ($gsi) - Changing transcript name from ${t_name}-$t_version to $new_name\n");
				my $xref_id = $trans_dbentry->dbID;
				my $num = $sth->execute($new_name, $xref_id) unless ($support->param('dry_run'));
			}
		}
	}
}


# finish logfile
$support->finish_log;

