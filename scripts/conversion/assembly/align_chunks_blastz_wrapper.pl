#!/usr/local/ensembl/bin/perl

=head1 NAME

align_chunks_blastz_wrapper.pl - lsf wrapper for align_chunks_blastz.pl

=head1 SYNOPSIS

align_chunks_blastz.pl [options]

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
    --evegadbname=NAME                  use ensembl-vega (target) database NAME
    --evegahost=HOST                    use ensembl-vega (target) database host
                                        HOST
    --evegaport=PORT                    use ensembl-vega (target) database port
                                        PORT
    --evegauser=USER                    use ensembl-vega (target) database
                                        username USER
    --evegapass=PASS                    use ensembl-vega (target) database
                                        passwort PASS
    --bindir=DIR                        look for program binaries in DIR

=head1 DESCRIPTION

This script runs align_chunks_blastz.pl one chromosome at a time via lsf on the
farm. See documention in align_chunks_blastz.pl for more information.

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
    'chromosomes|chr=s@',
    'evegahost=s',
    'evegaport=s',
    'evegauser=s',
    'evegapass=s',
    'evegadbname=s',
    'bindir=s',
);
$support->allowed_params(
    $support->get_common_params,
    'chromosomes',
    'evegahost',
    'evegaport',
    'evegauser',
    'evegapass',
    'evegadbname',
    'bindir',
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

# connect to database and get adaptors
my $E_dba = $support->get_database('evega', 'evega');
my $E_dbh = $E_dba->dbc->db_handle;

# get Vega chromosomes
my $V_chrlength = $support->get_chrlength($E_dba, $support->param('assembly'));

# loop over chromosomes
$support->log_stamped("Looping over chromosomes...\n");

foreach my $V_chr ($support->sort_chromosomes($V_chrlength)) {
    $support->log_stamped("Chromosome $V_chr...\n", 1);

    # run align_chunks_blastz.pl via lsf
    my $options = $support->create_commandline_options({
        'allowed_params' => 1,
        'replace' => {
            logfile     => "align_chunks_blastz.$V_chr.log",
            chromosomes => $V_chr,
            interactive => 0,
            bindir      => '/usr/local/ensembl/bin',
        },
    });
    $support->log("Running align_chunks_blastz.pl via lsf...\n", 2);

    my $logpath = $support->param('logpath').'/lsf';
    unless (-e $logpath) {
        system("mkdir $logpath") == 0 or
            $support->log_error("Can't create lsf log dir $logpath: $!\n");
    }
    my $cmd = qq(bsub -R'type=LINUX86' -o $logpath/align_chunks_blastz.%J.out -e $logpath/align_chunks_blastz.%J.err ./align_chunks_blastz.pl $options
    );

    system("$cmd") == 0
        or $support->throw("Error running align_chunks_blastz.pl: $!");

    $support->log_stamped("Done.\n", 2);

    $support->log_stamped("Done with chromosome $V_chr.\n", 1);
}
$support->log_stamped("Done.\n");

# finish logfile
$support->finish_log;

