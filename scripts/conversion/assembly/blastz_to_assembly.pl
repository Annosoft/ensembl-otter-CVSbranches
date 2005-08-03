#!/usr/local/ensembl/bin/perl

=head1 NAME

blastz_to_assembly.pl - create assembly from blastz whole genome alignment
between two closely related assemblies

=head1 SYNOPSIS

blastz_to_assembly.pl [options]

General options:
    --conffile, --conf=FILE             read parameters from FILE
                                        (default: conf/Conversion.ini)

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
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use File::Path;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'evegahost=s',
    'evegaport=s',
    'evegauser=s',
    'evegapass=s',
    'evegadbname=s',
);
$support->allowed_params(
    $support->get_common_params,
    'evegahost',
    'evegaport',
    'evegauser',
    'evegapass',
    'evegadbname',
);


if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
my $id = $ENV{'LSB_JOBINDEX'};
$support->param('logfile', "blastz_to_assembly.block${id}.log");
$support->init_log;

# log LSF execution host (so that we can track down errors in /tmp files)
$support->log("Running on host $ENV{LSB_HOSTS}.\n");

# connect to database and get adaptors
my ($V_dba, $V_dbh, $E_dba, $E_dbh);
$V_dba = $support->get_database('core');
$V_dbh = $V_dba->dbc->db_handle;
my $V_sa = $V_dba->get_SliceAdaptor;
$E_dba = $support->get_database('evega', 'evega');
$E_dbh = $E_dba->dbc->db_handle;
my $E_sa = $E_dba->get_SliceAdaptor;

# get non-aligned regions from tmp_align table
my $sth = $E_dbh->prepare(qq(SELECT * FROM tmp_align WHERE tmp_align_id = $id));
$sth->execute;
my $row = $sth->fetchrow_hashref;
$sth->finish;

# get slices from both assemblies
my $E_slice = $E_sa->fetch_by_region(
    'chromosome',
    $row->{'seq_region_name'},
    $row->{'e_start'},
    $row->{'e_end'},
    1,
    $support->param('ensemblassembly'),
);
my $V_slice = $V_sa->fetch_by_region(
    'chromosome',
    $row->{'seq_region_name'},
    $row->{'v_start'},
    $row->{'v_end'},
    1,
    $support->param('assembly'),
);

# write sequences to file
$support->log_stamped("Writing sequences to fasta and nib files...\n");
my $tmpdir = '/tmp/pm2.'.time;
system("mkdir $tmpdir") == 0 or
    $support->log_error("Can't create tmp dir $tmpdir: $!\n");
my $E_fh = $support->filehandle('>', "$tmpdir/e_seq.$id.fa");
my $V_fh = $support->filehandle('>', "$tmpdir/v_seq.$id.fa");
print $E_fh join(':', ">e_seq.$id dna:chromfrag chromosome",
                      $support->param('ensemblassembly'),
                      $E_slice->start,
                      $E_slice->end,
                      $E_slice->strand
                ), "\n";
print $E_fh $E_slice->get_repeatmasked_seq(undef, 1)->seq, "\n";
close($E_fh);
print $V_fh join(':', ">v_seq.$id dna:chromfrag chromosome",
                      $support->param('assembly'),
                      $V_slice->start,
                      $V_slice->end,
                      $V_slice->strand
                ), "\n";
print $V_fh $V_slice->get_repeatmasked_seq(undef, 1)->seq, "\n";
close($V_fh);

# convert sequence files from fasta to nib format (needed for lavToAxt)
system("faToNib $tmpdir/e_seq.$id.fa $tmpdir/e_seq.$id.nib") == 0 or
    $support->log_error("Can't run faToNib: $!\n");
system("faToNib $tmpdir/v_seq.$id.fa $tmpdir/v_seq.$id.nib") == 0 or
    $support->log_error("Can't run faToNib: $!\n");
$support->log_stamped("Done.\n");

# align using blastz
$support->log_stamped("Running blastz...\n");
my $blastz_cmd = qq(blastz $tmpdir/e_seq.$id.fa $tmpdir/v_seq.$id.fa Q=blastz_matrix.txt T=0 L=10000 H=2200 Y=3400 > $tmpdir/blastz.$id.lav);
system($blastz_cmd) == 0 or $support->log_error("Can't run blastz: $!\n");
$support->log_stamped("Done.\n");

# convert blastz output from lav to axt format
$support->log_stamped("Converting blastz output from lav to axt format...\n");
system("lavToAxt $tmpdir/blastz.$id.lav $tmpdir $tmpdir $tmpdir/blastz.$id.axt") == 0 or $support->log_error("Can't run lavToAxt: $!\n");
$support->log_stamped("Done.\n");

# find best alignment with axtBest
$support->log_stamped("Finding best alignment with axtBest...\n");
system("axtBest $tmpdir/blastz.$id.axt all $tmpdir/blastz.$id.best.axt") == 0 or $support->log_error("Can't run axtBest: $!\n");
$support->log_stamped("Done.\n");

# parse blastz output (see Graham's parser for chimp)

# write alignment to assembly table

# log alignment stats
# things to log:
#   - alignment length
#   - gap length
#   - number of aligned blocks

# cleanup
$support->log_stamped("Cleanup: Removing tmpdir...\n");
#rmtree($tmpdir) or $support->log_warning("Could not delete $tmpdir: $!\n");
$support->log_stamped("Done.\n");

# finish logfile
$support->finish_log;

