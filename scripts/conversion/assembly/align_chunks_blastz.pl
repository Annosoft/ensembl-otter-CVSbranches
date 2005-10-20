#!/usr/local/bin/perl

=head1 NAME

align_chunks_blastz.pl - create whole genome alignment between chunks of two
closely related assemblies using blastz from scratch

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
assemblies using blastz. It was originally written for zebrafish. Alignments
are created between 30Mb chunks of Ensembl assembly and max. 10Mb chunks of
annotated clones from the Vega assembly.

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

# log LSF execution host (so that we can track down errors in /tmp files)
if ($ENV{'LSB_HOSTS'}) {
    $support->log("Running on host $ENV{LSB_HOSTS}.\n");
    $support->log("LSF ID is $ENV{LSB_JOBID}.\n");
}

# connect to database and get adaptors
my $V_dba = $support->get_database('core');
my $V_dbh = $V_dba->dbc->db_handle;
my $V_sa = $V_dba->get_SliceAdaptor;
my $E_dba = $support->get_database('evega', 'evega');
my $E_dbh = $E_dba->dbc->db_handle;
my $E_sa = $E_dba->get_SliceAdaptor;

# create BlastzAligner object
my $aligner = BlastzAligner->new(-SUPPORT => $support);

# create tmpdir to store input and output
$aligner->create_tempdir;

# get Vega and Ensembl chromosomes
my $V_chrlength = $support->get_chrlength($E_dba, $support->param('assembly'));
my $E_chrlength = $support->get_chrlength($E_dba, $support->param('ensemblassembly'));
my $ensembl_chr_map = $support->get_ensembl_chr_mapping($V_dba, $support->param('assembly'));

# loop over chromosomes
$support->log_stamped("Looping over chromosomes...\n");

foreach my $V_chr ($support->sort_chromosomes($V_chrlength)) {
    $support->log_stamped("Chromosome $V_chr...\n", 1);
    $aligner->seq_region_name($V_chr);
    
    # skip non-ensembl chromosomes (e.g. MHC haplotypes)
    my $E_chr = $ensembl_chr_map->{$V_chr};
    unless ($E_chrlength->{$E_chr}) {
        $support->log("No equivalent chromosome in Ensembl. Skipping.\n", 1);
        next;
    }

    # get 30Mb chunks of Ensembl chromsomes
    my @E_slices;
    my $E_chr_slice = $E_sa->fetch_by_region('chromosome', $E_chr, undef,
        undef, undef, $support->param('ensemblassembly'));
    my $chunk_size = 30000000;
    my $start = 1;
    my $chunks = int($E_chr_slice->length/$chunk_size);
    for (my $i = 1; $i <= $chunks; $i++) {
        push @E_slices, $E_sa->fetch_by_region(
            'chromosome',
            $E_chr,
            $start,
            $i*$chunk_size,
            1,
            $support->param('ensemblassembly')
        );
        $start = $i*$chunk_size+1;
    }
    push @E_slices, $E_sa->fetch_by_region(
        'chromosome',
        $E_chr,
        $start,
        $E_chr_slice->length,
        1,
        $support->param('ensemblassembly')
    ) if ($start < $E_chr_slice->length);

    # get annotated chunks of Vega chromosomes
    my @V_slices;
    my $V_chr_slice = $V_sa->fetch_by_region('chromosome', $V_chr, undef,
        undef, undef, $support->param('assembly'));
    my @noannot = @{ $V_chr_slice->get_all_MiscFeatures('NoAnnotation') || [] };
    my @annot = &invert_selection($V_chr_slice->length, \@noannot);
    if (scalar(@annot)) {
        
        foreach my $coord (@annot) {
            push @V_slices, $V_sa->fetch_by_region(
                'chromosome',
                $V_chr,
                $coord->[0],
                $coord->[1],
                1,
                $support->param('assembly')
            );
        }
    } else {
        $support->log("No annotated regions on this chromosome. Skipping.\n", 1);
        next;
    }
    
    if ($support->param('verbose')) {
        my $fmt1 = "%-13s%-13s\n";
        my $fmt2 = "%11.0f  %11.0f\n";

        $support->log("\nEnsembl blocks:\n", 2);
        $support->log(sprintf($fmt1, qw(START END)), 3);
        $support->log(('-'x26)."\n", 3);
        foreach my $slice (@E_slices) {
            $support->log(sprintf($fmt2, $slice->start, $slice->end), 3);
        }
        
        $support->log("\nVega blocks:\n", 2);
        $support->log(sprintf($fmt1, qw(START END)), 3);
        $support->log(('-'x30)."\n", 3);
        foreach my $slice (@V_slices) {
            $support->log(sprintf($fmt2, $slice->start, $slice->end), 3);
        }

        $support->log("\n");
    }

    # write Vega sequences to file, and convert sequence files from fasta to
    # nib format (needed for lavToAxt)
    $support->log("Writing Vega sequences to fasta and nib files...\n", 2);
    my $V_basename = "v_seq.$V_chr";
    my %V_coords;
    my $i;
    foreach my $V_slice (@V_slices) {
        $i++;
        $aligner->write_sequence(
            $V_slice,
            $support->param('assembly'), 
            "$V_basename.$i",
            $V_basename
        );

        $V_coords{"$V_chr.$i"} = [$V_slice->start, $V_slice->end];
    }
    $support->log("Done.\n", 2);

    # loop over Ensembl blocks
    $support->log("Looping over Ensembl blocks...\n", 2);
    my $j;
    foreach my $E_slice (@E_slices) {
        $j++;
        my $E_basename = "e_seq.$V_chr.$j";
        $aligner->id("$V_chr.$j");

        # write Ensembl sequences to file, and convert sequence files from
        # fasta to nib format (needed for lavToAxt)
        $support->log("Writing Ensembl sequence to fasta and nib files...\n", 3);
        $aligner->write_sequence(
            $E_slice,
            $support->param('ensemblassembly'),
            $E_basename
        );
        $support->log("Done.\n", 3);

        # align using blastz
        $support->log("Running blastz...\n", 3);
        $aligner->run_blastz($E_basename, $V_basename);
        $support->log("Done.\n", 3);

        # convert blastz output from lav to axt format
        $support->log("Converting blastz output from lav to axt format...\n", 3);
        $aligner->lav_to_axt;
        $support->log("Done.\n", 3);

        # find best alignment with axtBest
        $support->log("Finding best alignment with axtBest...\n", 3);
        $aligner->find_best_alignment;
        $support->log("Done.\n", 3);

        # parse blastz output, and convert relative alignment coordinates to
        # chromosomal coords
        $support->log("Parsing blastz output...\n", 3);
        $aligner->parse_blastz_output;
        $aligner->adjust_coords(
            $E_slice->start,
            $E_slice->end,
            \%V_coords
        );
        $support->log("Done.\n", 3);

        # log alignment stats
        $aligner->log_block_stats(3);
    }
    $support->log("Done.\n", 2);

    $support->log_stamped("Done with chromosome $V_chr.\n", 1);
}
$support->log_stamped("Done.\n\n");

# filter overlapping Vega alignment regions
$support->log_stamped("Filtering overlapping Vega alignment regions...\n");
$aligner->filter_overlaps;
$support->log_stamped("Done.\n");

unless ($support->param('dry_run')) {
    # write alignments to assembly table
    $aligner->write_assembly($V_dba, $E_dbh, $E_sa);
}

# cleanup
$support->log_stamped("\nRemoving tmpdir...\n");
$aligner->remove_tempdir;
$support->log_stamped("Done.\n");

# overall stats
$aligner->log_overall_stats;

# finish logfile
$support->finish_log;


### END main ###

sub invert_selection {
    my ($chr_length, $noannot) = @_;

    # create start/end pairs of annotated regions (this is the inversion of the
    # non-annotated regions stored as MiscFeatures in the database)
    my $start = 1;
    my @annot_temp;
    foreach my $mf (@{ $noannot }) {
        push @annot_temp, [$start, $mf->start];
        $start = $mf->end;
    }
    push @annot_temp, [$start, $chr_length];

    # filter zero-length pairs, and split regions longer than 10Mb
    my @annot;
    my $chunk_size = 10000000;
    foreach my $coord (@annot_temp) {
        my $length = $coord->[1] - $coord->[0];
        next unless ($length > 0);

        my $start = $coord->[0];
        my $chunks = int($length/$chunk_size);
        for (my $i = 1; $i <= $chunks; $i++) {
            push @annot, [$start, $start+$i*$chunk_size];
            $start += $i*$chunk_size+1;
        }
        push @annot, [$start, $coord->[1]] if ($start < $coord->[1]);
    }

    return @annot;
}

