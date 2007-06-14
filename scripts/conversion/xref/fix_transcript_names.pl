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
    --prune=0|1                         remove changes from previous runs of this script
    -v, --verbose                       verbose logging (default: false)
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)

Specific options:

    --fix_xrefs, --fixxrefs=0|1         also fix xrefs

=head1 DESCRIPTION

This script updates transcript_info.name (which originally holds clone-name
based transcript names) with names that reflect the gene name.

Where the above results in identical names for transcripts from the same gene then
the transcripts are numbered incrementaly after ordering from the longest coding
to the shortest non-coding. The exceptions to this is are genes that have a
'Annotation_remark - fragmented loci' gene_remark, or a '%fragmented%' transcript_remark
 - these are truly fragmented genes so a remark to this effect is added and the transcript IDs are not patched. The list
of gene IDs at the end of the script are those that have been checked by Havana -
new ones are reported and should be passed back to Havana to check that they're not
true fragmented loci.

If you run the script after you ran add_vega_xrefs.pl, you can use the
--fix_xrefs option to also patch the diplay_xrefs. Note that if this option is used,
as it is might be with zebrafish, then transcript_info.name as well as the display_xref
will be updated to any imported Zfin name. Take care with this legacy option!

CAUTION: The -prune option restores the data to the stage before this script was *first*
run, ie do not use this option if add_vega_xrefs.pl has been run, or if gene or transcript_attribs
have been added between times. Or else be prepared to regenerate them.

It cannot be used to change the root of transcript names after the addition of Zfin gene
names - that must be done using patch_transcript_names.pl.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>
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
$support->parse_extra_options('fix_xrefs=s','prune=s');
$support->allowed_params($support->get_common_params, 'fix_xrefs','prune');

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

#get list of IDs that have previously been sent to annotators
my %seen_genes;
while (<DATA>) {
	next if /^\s+$/ or /#+/;
	$seen_genes{$_} = 1;
}

# connect to database and get adaptors
my $dba = $support->get_database('vega');
my $ga  = $dba->get_GeneAdaptor;
my $ta  = $dba->get_TranscriptAdaptor;
my $aa  = $dba->get_AttributeAdaptor;
my $dbh = $dba->dbc->db_handle;

# make a backup table the first time the script is run
my (%tabs);
map { $_ =~ s/`//g; $tabs{$_} += 1; } $dbh->tables;
if (! exists ($tabs{'backup_transcript_info_original'})) {
	$support->log("Creating backup of transcript_info table\n\n");
	$dbh->do("CREATE table backup_transcript_info_original
              SELECT * FROM transcript_info");
	$support->log("Creating backup of transcript_attrib\n\n");
	$dbh->do("CREATE table backup_transcript_attrib_original
              SELECT * FROM transcript_attrib");
	$support->log("Creating backup of gene_attrib\n\n");
	$dbh->do("CREATE table backup_gene_attrib_original
              SELECT * FROM gene_attrib");
	$support->log("Creating backup of xref\n\n");
	$dbh->do("CREATE table backup_xref_original
              SELECT * FROM xref");
}

if ($support->param('prune') 
		&& (exists ($tabs{'backup_transcript_info_original'}))
			&& $support->user_proceed("\nDo you want to undo changes from previous runs of this script?")) {
    $support->log("Undoing changes from previous runs of this script...\n");
	$dbh->do("DELETE from transcript_info");
	$dbh->do("INSERT into transcript_info SELECT * FROM backup_transcript_info_original");
	$dbh->do("DELETE from transcript_attrib");
	$dbh->do("INSERT into transcript_attrib SELECT * FROM backup_transcript_attrib_original");
	$dbh->do("DELETE from xref");
	$dbh->do("INSERT into xref SELECT * FROM backup_xref_original");
	$dbh->do("DELETE from gene_attrib");
	$dbh->do("INSERT into gene_attrib SELECT * FROM backup_gene_attrib_original");
}

#are duplicate transcript names going to be fixed later on ?
my $fix_names = 0;
if ($support->user_proceed("Do you want to check and fix duplicated transcript names?") ) {
	$fix_names = 1;
}

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
$support->log("Updating transcripts...\n");
my %stats = ('transcript_info' => 0, 'patched_transcript_info' =>0, 'xref' => 0, 'patched_xref' => 0);
my %transnames;
while (my ($gene_name, $gsi, $tsi, $trans_name) = $sth->fetchrow_array) {
    if ($trans_name =~ /(\-\d+)$/) {
        my $new_name = "$gene_name$1";
		$transnames{$gsi}->{$new_name}++;
        next if ($new_name eq $trans_name);
		
		#store a transcript attrib for the old name
        my $attrib = [
			Bio::EnsEMBL::Attribute->new(
                -CODE => 'synonym',
                -NAME => 'Synonym',
                -DESCRIPTION => 'Synonymous names',
                -VALUE => $trans_name,
            )
		];
		my $trans = $ta->fetch_by_stable_id($tsi);
		my $t_dbID = $trans->dbID;
		unless ($support->param('dry_run')) {
			$aa->store_on_Transcript($t_dbID, $attrib);
		}
		$support->log_verbose("Stored transcript attrib for old name for transcript $tsi\n");
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

#check for duplicated names and fix as appropriate
if ($fix_names) {
	foreach my $gid (keys %transnames) {
		my $gene;
		my %transcripts = %{$transnames{$gid}};
		if ( grep { $transcripts{$_} > 1 } keys %transcripts ) {
			$gene = $ga->fetch_by_stable_id($gid);
			&update_names($gene,\%stats,$dbh);
		}
	}
	$support->log("Patched $stats{patched_transcript_info} transcript_infos, $stats{patched_xref} xrefs.\n");
}

# finish logfile
$support->finish_log;

sub update_names {
	my ($gene,$stats,$dbh) = @_;
	my $gid    = $gene->stable_id;
	my $g_dbID = $gene->dbID;
	#get gene name from xref (if xrefs are to be fixed) or gene_name
	my $gene_name = ($support->param('fix_xrefs')) ? $gene->display_xref->display_id || $gene->gene_info->name->name
                    : $gene->gene_info->name->name;
	my @gene_remarks = $gene->gene_info->remark;

#	#get transcript_remarks
	my %trans_remarks;
	foreach my $trans (@{$ta->fetch_all_by_Gene($gene)}) {
		my @remarks = $trans->transcript_info->remark;
		foreach my $remark (@remarks) {
			my $r = $remark->remark;
			$trans_remarks{$r}++;
		}
	}

	#add a gene remark for the website if it's a framented gene
	if ( (grep { $_->remark eq 'Annotation_remark- fragmented_loci'} @gene_remarks )
			 || (grep { $_ =~ /fragmen/} keys %trans_remarks ) ) {
		my $comment = 'This locus has been annotated as fragmented because either there is not enough evidence covering the whole locus to identify the exact exon structure of the transcript, or because the transcript spans a gap in  the assembly';
		if (grep { $_->remark eq $comment} @gene_remarks) {
			$support->log("Fragmented loci annotation remark for gene $gid already exists\n");
		}
		else {
			my $attrib = [
				Bio::EnsEMBL::Attribute->new(
					-CODE => 'remark',
					-NAME => 'Remark',
					-DESCRIPTION => 'Annotation remark',
					-VALUE => $comment,
				)
				];
			unless ( $support->param('dry_run') ) {
				$aa->store_on_Gene($g_dbID,$attrib);
			}			
			$support->log("Added fragmented loci annotation remark for gene $gid\n");
		}
		return;
	}
	
	#otherwise patch the transcript_name  if not OKeyd before (and xref if requested)
	elsif (! $seen_genes{$gid} ) {
		$support->log_warning("$gid has duplicated names and has no \'Annotation_remark- fragmented_loci\' on the gene or \'\%fragmen\%\' remark on any transcripts. Transcript names are being patched but this ID should be reported to Havana for double checking\n");

		my @trans = $gene->get_all_Transcripts();
		#seperate coding and non_coding transcripts
		my $coding_trans = [];
		my $noncoding_trans = [];
		foreach my $trans ( @{$gene->get_all_Transcripts()} ) {
			if ($trans->translate) {
				push @$coding_trans, $trans;
			}
			else {
				push @$noncoding_trans, $trans;
			}
		}
		my $c = 0;
		#sort transcripts coding > non-coding, then on length
		foreach my $array_ref ($coding_trans,$noncoding_trans) {
			foreach my $trans ( sort { $b->length <=> $a->length } @$array_ref ) {
				my $tsi = $trans->stable_id;
				$c++;
				my $ext = sprintf("%03d", $c);
				my $new_name = $gene_name.'-'.$ext;
				$support->log("Gene $gid ($gene_name): transcript $tsi will be patched to $new_name\n"); 
				unless ($support->param('dry_run')) {
					# update transcript_info.name
					$stats->{'patched_transcript_info'} += $dbh->do(qq(
                                 UPDATE  transcript_info
                                 SET     name = "$new_name"
                                 WHERE   transcript_stable_id = "$tsi"));
					# fix xrefs too if asked for
					if ($support->param('fix_xrefs')) {
						$stats->{'patched_xref'} += $dbh->do(qq(
                                      UPDATE  xref x, external_db ed
                                      SET     x.display_label = "$new_name"
                                      WHERE   dbprimary_acc = "$tsi"
                                      AND     x.external_db_id = ed.external_db_id
                                      AND     ed.db_name = "Vega_transcript"));
					}					
				}
			}
		}		
	}			
}


#add details of genes with duplicated names that have already been reported to Havana...
__DATA__

OTTMUSG00000005478
OTTMUSG00000001936
OTTMUSG00000003440
OTTMUSG00000017081
OTTMUSG00000016310
OTTMUSG00000011441
OTTMUSG00000012302
OTTMUSG00000013368
OTTMUSG00000013335
OTTMUSG00000015766
OTTMUSG00000016025
OTTMUSG00000001066
OTTMUSG00000016331
OTTMUSG00000006935
OTTMUSG00000011654
OTTMUSG00000001835
OTTMUSG00000007263
OTTMUSG00000000304
OTTMUSG00000009150
OTTMUSG00000008023
OTTMUSG00000017077
