#!/usr/local/bin/perl

=head1 NAME

assembly_tags.pl - transform assembly_tags into MiscFeatures

=head1 SYNOPSIS

    assembly_tags.pl [options]

    General options:
        --dbname, db_name=NAME              use database NAME
        --host, --dbhost, --db_host=HOST    use database host HOST
        --port, --dbport, --db_port=PORT    use database port PORT
        --user, --dbuser, --db_user=USER    use database username USER
        --pass, --dbpass, --db_pass=PASS    use database passwort PASS
        --driver, --dbdriver, --db_driver=DRIVER    use database driver DRIVER
        --conffile, --conf=FILE             read parameters from FILE
        --logfile, --log=FILE               log to FILE (default: *STDOUT)
        -i, --interactive                   run script interactively
                                            (default: true)
        -n, --dry_run, --dry                don't write results to database
        -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script tranforms assembly_tags into MiscFeatures with MiscSet
'AssemblyTag'. Tag type and info are stored as 'assembly_tag' Attributes.

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
use Data::Dumper;

BEGIN {
    $SERVERROOT = "$Bin/../../../..";
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::Attribute;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->log_filehandle('>>');
$support->log($support->init_log);

# get adaptors
my $dba = $support->get_database('ensembl');
my $sa = $dba->get_SliceAdaptor;
my $mfa = $dba->get_MiscFeatureAdaptor;
my $msa = $dba->get_MiscSetAdaptor;
my $dbh = $dba->dbc->db_handle;

# create a MiscSet
$support->log("Creating MiscSet 'AssemblyTag'...\n");
my $mset = Bio::EnsEMBL::MiscSet->new(
   -CODE => 'AssemblyTag',
   -NAME => 'Assembly tag',
   -DESCRIPTION => 'Assembly tag',
   -LONGEST_FEATURE => 1e10,
);
$msa->store($mset) unless $support->param('dry_run');
$support->log("Done.\n");

# read data from assembly_tag
$support->log("Reading assembly_tag and storing as MiscFeature...\n");
my $sth = $dbh->prepare("SELECT * FROM assembly_tag");
$sth->execute;
my $i;
while my $row ($sth->fetchrow_hashref) {
    # create a slice
    my $slice = $sa->fetch_by_seq_region_id($row->{'SEQ_REGION_ID'});
    
    # create MiscFeature and add MiscSet
    my $mfeat = Bio::EnsEMBL::MiscFeature->new(
            -START  => $row->{'SEQ_REGION_START'},
            -END    => $row->{'SEQ_REGION_END'},
            -STRAND => $row->{'SEQ_REGION_STRAND'},
            -SLICE  => $slice,
    );
    $mfeat->add_MiscSet($mset);

    # create Attribute and add to MiscFeature
    my $attrib = Bio::EnsEMBL::Attribute->new(
            -CODE   => 'assembly_tag',
            -NAME   => $row->{'TAG_TYPE'},
            -VALUE  => $row->{'TAG_INFO'},
    );
    $mfeat->add_Attribute($attrib);

    # store MiscFeature
    $mfa->store($misc_feature) unless $support->param('dry_run');

    # stats
    $i++;
}
$support->log("Done storing $i assembly tags.\n");

# drop now obsolete table assembly_tag
$support->log("Dropping table assembly_tag...\n");
$dbh->do("DROP TABLE assembly_tag");
$support->log("Done.\n");

# finish log
$support->log($support->finish_log);


