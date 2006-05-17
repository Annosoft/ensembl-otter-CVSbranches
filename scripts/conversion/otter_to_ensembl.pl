#!/usr/local/bin/perl

=head1 NAME

otter_to_ensembl.pl - store Otter-specific data in core Ensembl tables

=head1 SYNOPSIS

otter_to_ensembl.pl [options]

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

This script stores all Vega-specific information needed for web display in core
Ensembl schema tables by putting data from gene_info/transcript_info and
related tables into gene_attrib and transcript_attrib. These mappings are done
at the moment:

    source                                  target              attrib code
    -----------------------------------------------------------------------
    author.author_name                      gene_attrib         author
    gene_synonym.name                       gene_attrib         synonym
    gene_remark.remark                      gene_attrib         remark
    author.author_name                      transcript_attrib   author
    transcript_remark.remark                transcript_attrib   remark
    transcript_info.cds_start_not_found     transcript_attrib   remark
    transcript_info.cds_end_not_found       transcript_attrib   remark
    transcript_info.mRNA_start_not_found    transcript_attrib   remark
    transcript_info.mRNA_end_not_found      transcript_attrib   remark

After running this script, the webcode no longer requires to load ensembl-otter.


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
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::EnsEMBL::Attribute;

$| = 1;

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'chromosomes|chr=s@',
);
$support->allowed_params(
    $support->get_common_params,
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

# connect to database and get adaptors
my $dba = $support->get_database('vega');
my $dbh = $dba->dbc->db_handle;
my $sa = $dba->get_SliceAdaptor;
my $aa = $dba->get_AttributeAdaptor;
my $ga = $dba->get_GeneAdaptor;

foreach my $chr ($support->sort_chromosomes) {
    $support->log_stamped("> Chr $chr\n");

    my $slice = $sa->fetch_by_region('chromosome', $chr);
    # loop over genes
    foreach my $gene (@{ $ga->fetch_by_Slice($slice) }) {
        $support->log_verbose($gene->display_xref->display_id." (".$gene->stable_id.")...\n", 1);
        
        # create attributes for all Vega-specific annotation
        my $gene_attribs = [];
        
        # author
        push @{ $gene_attribs }, Bio::EnsEMBL::Attribute->new(
            -CODE => 'author',
            -NAME => 'Author',
            -DESCRIPTION => 'Group resonsible for Vega annotation',
            -VALUE => $gene->gene_info->author->name
        );
        push @{ $gene_attribs }, Bio::EnsEMBL::Attribute->new(
            -CODE => 'author_email',
            -NAME => 'Author email address',
            -DESCRIPTION => 'Author\'s email address',
            -VALUE => $gene->gene_info->author->email
        );

        # gene_synonym
        foreach my $syn ($gene->gene_info->synonym) {
            push @{ $gene_attribs }, Bio::EnsEMBL::Attribute->new(
                -CODE => 'synonym',
                -NAME => 'Synonym',
                -DESCRIPTION => 'Synonymous names',
                -VALUE => $syn->name
            );
        }

        # gene_remark
        foreach my $remark ($gene->gene_info->remark) {
            # save remarks as long as they are not just whitespace!
            push @{ $gene_attribs }, Bio::EnsEMBL::Attribute->new(
                -CODE => 'remark',
                -NAME => 'Remark',
                -DESCRIPTION => 'Annotation remark',
                -VALUE => $remark->remark
            ) if ($remark->remark =~ /\S/);
        }

        # store attributes
        $support->log_verbose("Storing ".scalar(@$gene_attribs)." gene attributes.\n", 2);
        unless ($support->param('dry_run')) {
            $aa->store_on_Gene($gene->dbID, $gene_attribs);
        }

        # loop over transcripts
        foreach my $transcript (@{ $gene->get_all_Transcripts }) {
            $support->log_verbose($transcript->display_xref->display_id." (".$transcript->stable_id.")...\n", 2);
            
            my $trans_attribs = [];

            # author
            push @{ $trans_attribs }, Bio::EnsEMBL::Attribute->new(
                -CODE => 'author',
                -NAME => 'Author',
                -DESCRIPTION => 'Group resonsible for Vega annotation',
                -VALUE => $transcript->transcript_info->author->name
            );
            push @{ $trans_attribs }, Bio::EnsEMBL::Attribute->new(
                -CODE => 'author_email',
                -NAME => 'Author email address',
                -DESCRIPTION => 'Author\'s email address',
                -VALUE => $transcript->transcript_info->author->email
            );

            # transcript_remark
            foreach my $remark ($transcript->transcript_info->remark) {
                # save remarks as long as they are not just whitespace!
                push @{ $trans_attribs }, Bio::EnsEMBL::Attribute->new(
                    -CODE => 'remark',
                    -NAME => 'Remark',
                    -DESCRIPTION => 'Annotation remark',
                    -VALUE => $remark->remark
                ) if ($remark->remark =~ /\S/);
            }

            # start/end_not_found tags: store as remarks
            push @{ $trans_attribs }, Bio::EnsEMBL::Attribute->new(
                -CODE => 'remark',
                -NAME => 'Remark',
                -DESCRIPTION => 'Annotation remark',
                -VALUE => 'CDS start not found'
            ) if ($transcript->transcript_info->cds_start_not_found);
            push @{ $trans_attribs }, Bio::EnsEMBL::Attribute->new(
                -CODE => 'remark',
                -NAME => 'Remark',
                -DESCRIPTION => 'Annotation remark',
                -VALUE => 'CDS end not found',
            ) if ($transcript->transcript_info->cds_end_not_found);
            push @{ $trans_attribs }, Bio::EnsEMBL::Attribute->new(
                -CODE => 'remark',
                -NAME => 'Remark',
                -DESCRIPTION => 'Annotation remark',
                -VALUE => 'mRNA start not found',
            ) if ($transcript->transcript_info->mRNA_start_not_found);
            push @{ $trans_attribs }, Bio::EnsEMBL::Attribute->new(
                -CODE => 'remark',
                -NAME => 'Remark',
                -DESCRIPTION => 'Annotation remark',
                -VALUE => 'mRNA end not found',
            ) if ($transcript->transcript_info->mRNA_end_not_found);

            # store attributes
            $support->log_verbose("Storing ".scalar(@$trans_attribs)." transcript attributes.\n", 3);
            unless ($support->param('dry_run')) {
                $aa->store_on_Transcript($transcript->dbID, $trans_attribs);
            }
        }
    }

    $support->log_stamped("Done with chr $chr.\n\n");
}

# finish logfile
$support->finish_log;

