#!/usr/local/bin/perl

=head1 NAME

align_nonident_regions_wrapper.pl - lsf wrapper for align_nonident_regions.pl

=head1 SYNOPSIS

align_nonident_regions_wrapper.pl [arguments]

Required arguments:

    --dbname, db_name=NAME              database name NAME
    --host, --dbhost, --db_host=HOST    database host HOST
    --port, --dbport, --db_port=PORT    database port PORT
    --user, --dbuser, --db_user=USER    database username USER
    --pass, --dbpass, --db_pass=PASS    database passwort PASS
    --assembly=ASSEMBLY                 assembly version ASSEMBLY

    --altdbname=NAME                    alternative database NAME
    --altassembly=ASSEMBLY              alternative assembly version ASSEMBLY

Optional arguments:

    --chromosomes, --chr=LIST           only process LIST chromosomes
    --bindir=DIR                        look for program binaries in DIR
    --tmpfir=DIR                        use DIR for temporary files (useful for
                                        re-runs after failure)

    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)

    --logfile, --log=FILE               log to FILE (default: *STDOUT)
    --logpath=PATH                      write logfile to PATH (default: .)
    --logappend, --log_append           append to logfile (default: truncate)

    -v, --verbose=0|1                   verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script is a wrapper around align_nonident_regions.pl to run one chromosome
at a time via lsf. See documentation in align_nonident_regions.pl for details.

=head1 RELATED SCRIPTS

The whole Ensembl-vega database production process is done by these scripts:

    ensembl-otter/scripts/conversion/assembly/make_ensembl_vega_db.pl
    ensembl-otter/scripts/conversion/assembly/align_by_clone_identity.pl
    ensembl-otter/scripts/conversion/assembly/align_nonident_regions.pl
    ensembl-otter/scripts/conversion/assembly/map_annotation.pl
    ensembl-otter/scripts/conversion/assembly/finish_ensembl_vega_db.pl

See documention in the respective script for more information.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <meidl@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

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
use BlastzAligner;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'assembly=s',
    'altdbname=s',
    'altassembly=s',
    'bindir=s',
    'tmpdir=s',
    'chromosomes|chr=s@',
);
$support->allowed_params(
    $support->get_common_params,
    'assembly',
    'altdbname',
    'altassembly',
    'bindir',
    'tmpdir',
    'chromosomes',
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

$support->check_required_params(
    'assembly',
    'altdbname',
    'altassembly'
);

#####
# connect to database and get adaptors
#
my $R_dba = $support->get_database('ensembl');

# loop over chromosomes
$support->log_stamped("Looping over chromosomes...\n");

foreach my $chr ($support->sort_chromosomes) {
    $support->log_stamped("Chromosome $chr...\n", 1);

    # run align_nonident_regions.pl via lsf
    my $options = $support->create_commandline_options({
        'allowed_params' => 1,
        'replace' => {
            logfile     => "align_nonident_regions.chr.$chr.log",
            chromosomes => $chr,
            interactive => 0,
        },
    });
    $support->log("Running align_nonident_regions.pl via lsf...\n", 2);

    my $logpath = $support->param('logpath').'/lsf';
    unless (-d $logpath) {
        system("mkdir $logpath") == 0 or
            $support->log_error("Can't create lsf log dir $logpath: $!\n");
    }

    my $cmd = qq(bsub -o $logpath/align_nonident_regions.$chr.%J.out -e $logpath/align_nonident_regions.$chr.%J.err ./align_nonident_regions.pl $options
    );

    system("$cmd") == 0
        or $support->log_error("Error running align_nonident_regions.pl: $!");

    $support->log_stamped("Done.\n", 2);

    $support->log_stamped("Done with chromosome $chr.\n", 1);
}

$support->log_stamped("Done.\n");

# finish logfile
$support->finish_log;


