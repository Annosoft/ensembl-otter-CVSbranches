#!/usr/local/bin/perl

=head1 NAME

fix_transcript_names.pl - update transcript names with external names from genes

=head1 SYNOPSIS

fix_transcript_names.pl [options]

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

    --fix_xrefs, --fixxrefs=0|1         also fix xrefs

=head1 DESCRIPTION

This script updates transcript_info.name (which originally holds clone-name
based transcript names) with names that reflect the gene name.

If you run the script after you ran add_vega_xrefs.pl, you can use the
--fix_xrefs option to also patch the diplay_xrefs. It cannot be used to change
transcript names after the addition of Zfin gene names - that must be done using
patch_transcript_names.pl

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
$support->parse_extra_options('fix_xrefs=s');
$support->allowed_params($support->get_common_params, 'fix_xrefs');

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
my $dbh = $dba->dbc->db_handle;

# get gene name to transcript_stable_id mapping from db
$support->log("Fetching gene name to transcript_stable_id mapping from db...\n");
my $sql = qq(
    SELECT  gn.name,
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
$support->log("Updating transcripts...\n");
my %stats = ('transcript_info' => 0, 'xref' => 0);
while (my ($gene_name, $tsi, $trans_name) = $sth->fetchrow_array) {
    if ($trans_name =~ /(\-\d+)$/) {
        my $new_name = "$gene_name$1";
        next if ($new_name eq $trans_name);

        $support->log_verbose(sprintf("%-20s%-3s%-20s", $trans_name, "->", $new_name)."\n", 1);

        # update transcript_info.name
        $stats{'transcript_info'} += $dbh->do(qq(
            UPDATE  transcript_info
            SET     name = "$new_name"
            WHERE   transcript_stable_id = "$tsi"
        )) unless ($support->param('dry_run'));

        # fix xrefs too if asked for
        if ($support->param('fix_xrefs')) {
            $stats{'xref'} += $dbh->do(qq(
                UPDATE  xref x, external_db ed
                SET     x.display_label = "$new_name"
                WHERE   dbprimary_acc = "$tsi"
                AND     x.external_db_id = ed.external_db_id
                AND     ed.db_name = "Vega_transcript"
            )) unless ($support->param('dry_run'));
        }
    } else {
        $support->log_warning("Not updating transcript name $trans_name: unexpected transcript name.\n", 1);
    }
}
$support->log("Done updating $stats{transcript_info} transcript_infos, $stats{xref} xrefs.\n");

# finish logfile
$support->finish_log;

