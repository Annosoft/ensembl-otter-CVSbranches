#!/usr/local/bin/perl

=head1 NAME

find_GD_Havana_loci.pl - find where a single loci has both GD and havana transcript names

=head1 SYNOPSIS

find_GD_Havana_loci.pl [options]

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
    --prune=0|1                         remove changes from previous runs of this script
    -v, --verbose                       verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)


=head1 DESCRIPTION

Quick script to look for loci that are a mixture of Havana and non-Havana transcripts (as
evidenced by the presence or abscence of a prefix such as GD:)

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
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::Otter::GeneRemark;
use Data::Dumper;

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
my $ga  = $dba->get_GeneAdaptor;
my $ta  = $dba->get_TranscriptAdaptor;
my $aa  = $dba->get_AttributeAdaptor;
my $dbh = $dba->dbc->db_handle;


# get gene name to transcript_stable_id mapping from db
$support->log("Fetching gene name to transcript_stable_id mapping from db...\n");
my $sql = qq(
    SELECT  gn.name,
            gi.gene_stable_id,
            ti.transcript_stable_id,
            ti.name 
    FROM    gene_name gn,
            gene_info gi,
            gene_stable_id gsi, 
            transcript t,
            transcript_stable_id tsi,
            transcript_info ti 
    WHERE   gn.gene_info_id = gi.gene_info_id 
    AND     gi.gene_stable_id = gsi.stable_id
    AND     gsi.gene_id = t.gene_id 
    AND     t.transcript_id = tsi.transcript_id 
    AND     tsi.stable_id = ti.transcript_stable_id
);
my $sth = $dbh->prepare($sql);
$sth->execute;
$support->log("Done.\n");

# loop over transcripts
my ($not_first,$tot_c,$prefixed_c,$noprefixed_c) = 0;
my $current_name;
my @names;
my ($gene_name, $gsi, $tsi, $trans_name, $last_gsi);
while (($gene_name, $gsi, $tsi, $trans_name) = $sth->fetchrow_array) {

	#ignore external genes
#    next if ($gene_name =~ /:/); #ignore all non Havana genes
    next if ($gene_name !~ /^GD:/);

	if ($gene_name eq $current_name) {
		push @names, $trans_name;
	}
	elsif ($not_first) {
		my ($tot_c,$prefixed_c);
		foreach my $name (@names) {
			$tot_c++;
			# see if the name has a prefix
			if ($name =~ /.+:.+/) {
				$prefixed_c++;
			}
			else {
				$noprefixed_c++;
			}	
		}
		my $names = join q{ }, @names;
		if ($prefixed_c && $noprefixed_c ) {
			$support->log_warning("Have a look at gene $last_gsi - transcript names are --$names--\n");
		}
		else {
			$support->log("Gene $last_gsi is OK - transcript names are --$names--\n");
		}

		#reset
		@names = ($trans_name);
		($tot_c,$prefixed_c,$noprefixed_c) = 0;
		$current_name = $gene_name;
		$last_gsi = $gsi;
	}
	else {
		push @names, $trans_name;
		$current_name = $gene_name;
		$last_gsi = $gsi;
		$not_first = 1;
	}

}			

foreach my $name (@names) {
	$tot_c++;
	# see if the name has a prefix
	if ($name =~ /.+:.+/) {
		$prefixed_c++;
	}
}
if ($prefixed_c && $noprefixed_c ) {
	my $names = join q{,}, @names;
	$support->log("Have a look at gene $last_gsi - transcript names are --$names--\n");
}

