#!/usr/local/bin/perl

=head1 NAME

align_chunks_blastz.pl - create whole genome alignment between chunks of two
closely related assemblies using blastz

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

This script is part of a series of scripts to transfer annotation from a
Vega to an Ensembl assembly. See "Related scripts" below for an overview of the
whole process.

It creates a whole genome alignment between chunks of two closely related
assemblies using blastz.

=head1 RELATED SCRIPTS

The whole Ensembl-vega database production process is done by these scripts:

    ensembl-otter/scripts/conversion/assembly/make_ensembl_vega_db.pl
    
    (
        ensembl-otter/scripts/conversion/assembly/align_by_clone_identity.pl
    and
        ensembl-otter/scripts/conversion/assembly/align_nonident_regions.pl
    )
    or
        ensembl-otter/scripts/conversion/assembly/align_chunks_blastz.pl

    ensembl-otter/scripts/conversion/assembly/map_annotation.pl
    ensembl-otter/scripts/conversion/assembly/finish_ensembl_vega_db.pl

See documention in the respective script for more information.

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

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

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
            binpath     => '/usr/local/ensembl/bin',
        },
    });
    $support->log("Running align_chunks_blastz.pl via lsf...\n", 2);

    my $cmd = qq(bsub -P vega
                    -R'select[LINUX86]' \
                    -o align_chunks_blastz.%J.out \
                    -e align_chunks_blastz.%J.err
                 align_chunks_blastz.pl $options
    );

    system($cmd) == 0
        or $support->throw("Error running align_chunks_blastz.pl: $!");

    $support->log_stamped("Done.\n", 2);

    $support->log_stamped("Done with chromosome $V_chr.\n", 1);
}
$support->log_stamped("Done.\n");

# finish logfile
$support->finish_log;

