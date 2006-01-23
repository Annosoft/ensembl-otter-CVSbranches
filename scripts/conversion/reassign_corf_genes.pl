#!/usr/local/bin/perl

=head1 NAME

reassign_corf_genes.pl - change analysis for CORF genes on finished clones

=head1 SYNOPSIS

reassign_corf_genes.pl [options]

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

=head1 DESCRIPTION

This script sets the analysis of CORF genes to 'otter' if they are on a clone
marked as annotated.

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
    $SERVERROOT = "$Bin/../../..";
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
$support->parse_extra_options('chromosomes|c=s@');
$support->allowed_params($support->get_common_params, 'chromosomes');

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# make sure annotation_status.pl was run before this script
exit unless $support->user_proceed("This script must run after annotation_status.pl. Have you run it?");

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# connect to database and get adaptors
my $dba = $support->get_database('ensembl');
my $sa = $dba->get_SliceAdaptor;
my $dbh = $dba->dbc->db_handle;

# SQL statement to change analysis
my $sql = qq(
    UPDATE gene g, analysis a
    SET g.analysis_id = a.analysis_id
    WHERE g.gene_id = ?
    AND a.logic_name = 'otter'
);
my $sth = $dbh->prepare($sql);

# loop over chromosomes
foreach my $chr ($support->sort_chromosomes) {
    $support->log_stamped("> Chr $chr\n\n");
    
    # get all CORF genes on this chromosome
    my $slice = $sa->fetch_by_region('chromosome', $chr);
    my @genes = @{ $slice->get_all_Genes('otter_corf') };
    my $i = 0;

    foreach my $gene (@genes) {
        # project gene onto clones; if any of them are annotated, reassign the
        # gene to analysis 'otter'
        my $is_annotated = 0;

        CLONE:
        foreach my $segment (@{ $gene->feature_Slice->project('clone') || [] }) {
            my $noannot = 0;
            my $cl_slice  = $segment->to_Slice;
            foreach my $at (@{ $cl_slice->get_all_Attributes('NoAnnotation') }) {
                $noannot = 1 if ($at->value);
            }

            unless ($noannot) {
                $is_annotated = 1;
                last CLONE;
            }
        }

        if ($is_annotated) {
            $support->log("Reassigning ".$gene->display_xref->display_id. " (".$gene->stable_id.", ".$gene->dbID.").\n", 1);
            $sth->execute($gene->dbID) unless ($support->param('dry_run'));
            $i++;
        }
    }
    
    $support->log("\nReassigned $i of ".scalar(@genes)." genes.\n");
    $support->log_stamped("Done with chromosome $chr.\n\n");
}

# finish logfile
$support->finish_log;

