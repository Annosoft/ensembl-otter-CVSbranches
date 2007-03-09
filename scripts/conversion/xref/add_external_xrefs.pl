#!/usr/local/bin/perl

=head1 NAME

add_external_xrefs.pl - adds xrefs to external databases from various types
of input files

=head1 SYNOPSIS

add_external_xrefs.pl [options]

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
    --gene_stable_id, --gsi=LIST|FILE   only process LIST gene_stable_ids
                                        (or read list from FILE)
    --xrefformat=FORMAT                 input file format FORMAT
                                        (hugo|tcag|imgt|refseq)
    --hugofile=FILE                     read Hugo input from FILE
    --tcagfile=FILE                     read TCAG input from FILE
    --imgtfile=FILE                     read IMGT input from FILE
    --refseqfile=FILE                   read Refseq input from FILE
    --ensembl_xref                      read xrefs from an Ensembl db
    --mismatch                          correct case mismatches in the db
                                          overrides dry-run, doesn't add xrefs
    --prune                             reset to the state before running this
                                        script (i.e. after running
                                        add_vega_xrefs.pl)

=head1 DESCRIPTION

This script parses input files from various sources - HUGO, MGI, TCAG (human chr 7 annotation),
IMGT (human major histocompatibility complex nomenclature), RefSeq and an Ensembl
database - and adds xrefs to the databases covered by the respective input source. If
appropriate the display names of genes are set accordingly.

It's worthwhile running the script first with -dry_run and -mismatch options to fix any
case errors in the Vega Gene names. Then run it normally to add xrefs. Note that if any
gene_names are found to have case errors then the transcript names must also be updated
using patch_transcript_names.pl.

HUGO and an Ensembl db can be used at the same time for human - this might be expanded in
the future to use locuslink files. TCAG (The Centre for Applied Genomics) is used exclusively
for human to add xrefs for externally annotated (Sick-Kids) genes on human chr 7. IMGT is
used to add xrefs for HLA genes on human haplotypes

For mouse two files are parsed - the first (mgivega) associates OTTMUS IDs with MGI marker 
symbols, the second (mgi) adds links to external databases.

Currently, these input formats are supported:

    hugo        => http://www.gene.ucl.ac.uk/nomenclature/data/gdlw_index.html
                   ('All data' in text format)
    mgivega     => ftp://ftp.informatics.jax.org/pub/reports/MGI_VEGA.rpt
    mgi         => ftp://ftp.informatics.jax.org/pub/reports/MGI_MouseHumanSequence.rpt
    tcag        => http://www.chr7.org/download/Dec2005/TCAG_ANNOTATION.gff.gz
    imgt        => by email Steven Marsh <marsh@ebi.ac.uk>
    ensemblxref => use core ensembl database

There are parsers for locuslink and refseq but these have not been tested recently
    refseq      => ftp://ftp.ncbi.nih.gov/genomes/__SPECIES__/RNA/rna.gkb

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>
Original code by Patrick Meidl <meidl@ebi.ac.uk> and Tim Hubbard <th@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use FindBin qw($Bin);
use vars qw($SERVERROOT);
use Storable;

BEGIN {
    $SERVERROOT = "$Bin/../../../..";
    unshift(@INC, "$SERVERROOT/ensembl/modules");
    unshift(@INC, "$SERVERROOT/bioperl-live");
}

use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::ConversionSupport;
use Bio::SeqIO::genbank;
use Data::Dumper;

$| = 1;

our $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'chromosomes|chr=s@',
    'gene_stable_id|gsi=s@',
    'xrefformat=s@',
    'hugofile=s',
    'mgivegafile=s',
    'mgifile=s',
    'refseqfile=s',
	'tcagfile=s',
	'imgtfile=s',
    'ensembl_xref',
    'ensemblhost=s',
    'ensemblport=s',
    'ensembluser=s',
    'ensemblpass=s',
    'ensembldbname=s',
    'mismatch',
    'prune',
);
$support->allowed_params(
    $support->get_common_params,
    'chromosomes',
    'gene_stable_id',
    'xrefformat',
    'hugofile',
    'mgivegafile',
    'mgifile',
    'refseqfile',
    'tcagfile',
    'imgtfile',
    'ensembl_xref',
    'ensemblhost',
    'ensemblport',
    'ensembluser',
    'ensemblpass',
    'ensembldbname',
    'mismatch',
    'prune',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');
#$support->comma_to_list('xrefformat');
$support->list_or_file('gene_stable_id');

# ask user to confirm parameters to proceed
$support->confirm_params;

# make sure the corrct sequence of options have been run
if ( $support->param('xrefformat') eq 'mgi' ) {
	exit unless $support->user_proceed("MGI files can be parsed only after first using the mgivega option. Have you done this?");
}
else {
	exit unless $support->user_proceed("This script must run after add_vega_xrefs.pl. Have you run it?");
}

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# check that --mismatch is combined with --dry_run
if ($support->param('mismatch')) {
    $support->log("--mismatch is set, therefore setting --dry_run to 1...\n");
    $support->param('dry_run', 1);
}

# connect to database and get adaptors (caching features on one slice only)
# get an ensembl database for better performance (no otter tables are needed)
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;
my $sa = $dba->get_SliceAdaptor();
my $ga = $dba->get_GeneAdaptor();
my $ea = $dba->get_DBEntryAdaptor();

# statement handle for display_xref_id update
my $sth_display_xref = $dba->dbc->prepare("UPDATE gene SET display_xref_id=? WHERE gene_id=?");

# statement handles for fixing case errors
my $sth_case1 = $dba->dbc->prepare("UPDATE gene_name set name = ? WHERE name = ?");
my $sth_case2 = $dba->dbc->prepare("UPDATE xref set display_label = ? WHERE display_label = ?");

#make sure mgivega xrefs aren't pruned by mistake
if ( ($support->param('prune')) && ($support->param('xrefformat') eq 'mgi') ) {
	if ($support->user_proceed('You should not be pruning mgi_vega xrefs now. Shall I switch off pruning for you?') ) {
		$support->param('prune',0);
	}
	else {
		exit;
	}
}

# delete all external xrefs if --prune option is used; basically this resets xrefs to the state before running this script
if ($support->param('prune') and $support->user_proceed('Would you really like to delete all *external* xrefs before running this script?')) {
	my $num;
	# xrefs
	$support->log("Deleting all external xrefs...\n");
	$num = $dba->dbc->do(qq(
           DELETE x
           FROM xref x, external_db ed
           WHERE x.external_db_id = ed.external_db_id
           AND ed.db_name not in ('Vega_gene', 'Vega_transcript',
                                  'Vega_translation', 'Interpro', 'CCDS', 'ENST')
        ));
	$support->log("Done deleting $num entries.\n");

	# object_xrefs
	$support->log("Deleting orphan object_xrefs...\n");
	$num = $dba->dbc->do(qq(
           DELETE ox
           FROM object_xref ox
           LEFT JOIN xref x ON ox.xref_id = x.xref_id
           WHERE x.xref_id IS NULL
        ));
	# gene.display_xref_id
	$support->log("Resetting gene.display_xref_id...\n");
	$num = $dba->dbc->do(qq(
           UPDATE gene g, gene_stable_id gsi, xref x
           SET g.display_xref_id = x.xref_id
           WHERE g.gene_id = gsi.gene_id
           AND gsi.stable_id = x.dbprimary_acc
        ));

	$support->log("Done deleting $num entries.\n");
}

my @gene_stable_ids = $support->param('gene_stable_id');
my %gene_stable_ids = map { $_, 1 } @gene_stable_ids;
my $chr_length = $support->get_chrlength($dba);
my @chr_sorted = $support->sort_chromosomes($chr_length);

# get species name
my $species = $support->get_species_scientific_name($dba);
my %species_lookup = ( 'Mus musculus' => 'Mm', 'Homo sapiens' => 'Hs' );
my $sp = $species_lookup{$species};

# parse input file
no strict 'refs';
my %primary;
my $ens_xrefs = {};
my $xrefs = {};
my @xref_sources;
my $lcmap = {};
my $xref_file = $SERVERROOT.'/'.$support->param('dbname').'-parsed_records.file';
my $lc_xref_file = $SERVERROOT.'/'.$support->param('dbname').'-lc-parsed_records.file';

#only look at Ensembl db if there is no input file to be parsed, or if the file to be parsed is hugo
if (
     (!$support->param('xrefformat'))
     || ($support->param('xrefformat') eq 'hugo')
   ) {
	if ($support->param( 'ensembl_xref' )) {
		if (-e $xref_file) {
			if ($support->user_proceed("Read E! xref records from a previously saved files - $xref_file ?\n")) {
				$ens_xrefs = retrieve($xref_file);
				$lcmap = retrieve($lc_xref_file);
			}
			else {
				$support->log_stamped("Ignoring existing input file and re-reading from E! database...\n");
				&parse_ensdb($ens_xrefs,$lcmap,$sp);
				$support->log_stamped("Finished reading from E! database, storing to files...\n");
				store($ens_xrefs,$xref_file);
				store($lcmap,$lc_xref_file);
			}
		}
		else {
			$support->log_stamped("No input files read, using E! database\n");
			&parse_ensdb($ens_xrefs,$lcmap,$sp);
			$support->log_stamped("Finished reading from E! database, storing to files...\n");
			store($ens_xrefs,$xref_file);			
			store($lcmap,$lc_xref_file);
		}
		push @xref_sources, [ 'ensembl', $ens_xrefs ];
	}
}

#parse the input file (independently of E! db parsing)
if ($support->param('xrefformat')){
	$support->log_stamped("Reading xref input files...\n");
	my $format = ($support->param('xrefformat'));
	my $parser = "parse_$format";
	&$parser($xrefs, $lcmap);
	push @xref_sources, [ $format, $xrefs ];
}

#warn Dumper(\@xref_sources);
#warn Dumper($lcmap);
#exit;

use strict 'refs';
$support->log_stamped("Done.\n\n");

# define each type of xref that can be set, and whether to set as display_xref or not
my %extdb_def = (
    HUGO                => ['KNOWN'    , 1],
    EntrezGene          => ['KNOWNXREF', 0],
    MarkerSymbol        => ['KNOWNXREF', 1],
    RefSeq_dna          => ['KNOWN'    , 0],
    RefSeq_peptide      => ['KNOWN'    , 0],
    MIM_GENE            => ['KNOWNXREF', 0],
    'Uniprot/SWISSPROT' => ['KNOWN'    , 0],
	TCAG                => ['XREF'     , 0],
	IMGT                => ['XREF'     , 0],
    Ens_Mm_gene         => ['XREF'     , 0],
	Ens_Hs_gene         => ['XREF'     , 0],
);

# loop over chromosomes
$support->log("Looping over chromosomes: @chr_sorted\n\n");
my $seen_xrefs;
my %overall_stats;
foreach my $chr (@chr_sorted) {
    $support->log_stamped("> Chromosome $chr (".$chr_length->{$chr}."bp).\n\n");

    # fetch genes from db
    $support->log("Fetching genes...\n");
    my $slice = $sa->fetch_by_region('chromosome', $chr);
    my $genes = $ga->fetch_all_by_Slice($slice);
    $support->log_stamped("Done fetching ".scalar @$genes." genes.\n\n");

    # loop over genes
    my %stats = map { $_ => 0 } keys %extdb_def;
    my %warnings = (
        wrong_case      => 0,
        nomatch_clone   => 0,
        nomatch_missing => 0,
    );
    my $gnum = 0;
 GENE:
	foreach my $gene (@$genes) {
		my $gsi = $gene->stable_id;
        my $gid = $gene->dbID;

		# filter to user-specified gene_stable_ids
        if (scalar(@gene_stable_ids)){
            next unless $gene_stable_ids{$gsi};
        }
		
        # catch missing display_xrefs here!!
        my $disp_xref = $gene->display_xref;
        my ($stripped_name,$gene_name,$xref_dbname,$prefix);
        if ($disp_xref) {
            $gene_name = $disp_xref->display_id;
			$xref_dbname = $disp_xref->dbname;
        } else {
            $support->log_warning("No display_xref found for gene $gid ($gsi). Skipping.\n");
            next GENE;
        }

        $support->log("Gene $gene_name ($gid, $gsi)...\n");
	
        # see if the gene_name has a prefix
        if ( ($prefix,$stripped_name) = $gene_name  =~ /(.*?):(.*)/) {
			unless ($stripped_name) {
				$stripped_name = $prefix;
				$prefix = 0;
			}
		}
		else { $stripped_name = $gene_name; }


		#skip if we're adding CTAG xrefs and this is not an SK gene
		next if ( ($support->param('xrefformat') eq 'tcag') && $prefix ne 'SK');		

		$gnum++;

		# get all names including synonyms
		push my (@gene_names), $stripped_name;
		push my (@lc_names), lc($stripped_name);
		my @syn_names;
		if (my @syns = @{$gene->get_all_Attributes('synonym')} ) {
			foreach my $syn (@syns) {
				my ($pref,$name) = $syn->value =~ /(.*?):(.*)/;
				push @syn_names, $syn->value;
				push @gene_names, $name;
				push @lc_names, lc($name);
			}
		}

		#look only for stable_ids if we are adding IMGT xrefs or MarkerSymbol ones the first time around
		if ($support->param('xrefformat') =~ /imgt|mgivega/) {
			@gene_names = ( $gsi );
			@lc_names = ( );
		}

		#get a list of all db_names for xrefs on this gene
		my %existing_dbnames;
		my $xrefs = $gene->get_all_DBEntries;
		foreach my $xref (@{$xrefs}){
			my $dbname = $xref->dbname;
			$existing_dbnames{$dbname} = 1;
		}

		#use previously set MarkerSymbol xrefs as searchable names if we're setting Marker Symbol xrefs for the second time
		if ( $support->param('xrefformat') eq 'mgi') {
			foreach my $xref (@{$xrefs}){
				if ($xref->dbname eq 'MarkerSymbol') {
					my $mgi_name = $xref->display_id;
					unless (grep {$_ eq $mgi_name} @gene_names) {
						push @gene_names, $xref->display_id;
					}
				}
			}
		}

		my $xref_found = 0;
		foreach my $name ( @gene_names) {
			next if $xref_found;
			$support->log("Searching for name $name\n",1);

			foreach my $source (@xref_sources) {
				next if $xref_found;
				$support->log_verbose("Examining source $source->[0]\n",1); 
				# do we have a match on gene_name
				if ($source->[1]->{$name}) {
					$xref_found = 1;
					#report withdrawn records (HUGO)
					if ($source->[0] eq 'hugo') {
						if ($source->[1]->{$name}->{'Approved Symbol'} =~ /~withdrawn/) {
							$support->log_warning("Vega name: $name matches a withdrawn HGNC record\n");
						}
						if (my $alias_for = $source->[1]->{$name}->{'Is an alias for'}){
							$support->log_warning("Vega name: $name is an alias / previous record for HGNC: $alias_for\n");
						}
					}
				DB: foreach my $extdb (keys %extdb_def) {

						#don't go any further if this gene already has an xref for this source
						if ($existing_dbnames{$extdb}) {
							$support->log("$extdb xref previously set for gene $gid, not storing a new one.\n", 1);
							next DB;
						}

						#go through each xref for this source
						my ($xid,$pid,$concat_xid);
						foreach my $concat_xid ( @{$source->[1]->{$name}->{$extdb}} ) {
	
							#catch empty xrefs
							next DB if (!$concat_xid || $concat_xid =~ /^\|\|$/);

							#get display and accession idss
							if ($concat_xid =~ /\|\|/) {
								($xid,$pid) = split /\|\|/, $concat_xid;
							}
							else {
								$pid = $xid = $concat_xid;						
							}
							$stats{$extdb}++;
							if (!$support->param('dry_run')) {

								#use an existing xref if there is one...
								my ($existing_xref,$dbID);
								if ($existing_xref = $ea->fetch_by_db_accession($extdb,$pid)) {
									$support->log("Using previous xref for gene $gid ($extdb display_id $xid, pid = $pid).\n", 1);
									$gene->add_DBEntry($existing_xref);
									$dbID = $ea->store($existing_xref, $gid, 'gene');
								}

								#... or else create a new one
								else {
									my $dbentry = Bio::EnsEMBL::DBEntry->new(
																			 -primary_id => $pid,
																			 -display_id => $xid,
																			 -version    => 1,
																			 -release    => 1,
																			 -dbname     => $extdb,
																			);
									$dbentry->status($extdb_def{$extdb}->[0]);
									$gene->add_DBEntry($dbentry);
									$dbID = $ea->store($dbentry, $gid, 'gene');
								}
								if ($dbID) {
									$support->log("Stored $extdb display_id $xid, pid = $pid for gene $gid (dbID $dbID).\n", 1);

									#do we want to update the display_xref ?
									if ($extdb_def{$extdb}->[1]) {

										#hack to set MGI as a display xref only if the vega name is the same as MGI one
										if ($support->param('xrefformat') eq 'mgivega') {
											if (lc($stripped_name) eq lc($xid) ) {
												if ($stripped_name ne $xid) {
													$support->log_warning("Not setting a display_xref using MGI record ($xid) - different case from the Vega gene $gid($stripped_name)\n",1);
													next GENE;
												}
												else {
													$support->log("Setting a display_xref - Vega gene $gid($stripped_name) has the same name as the MGI record ($xid)\n",1);
#													$sth_display_xref->execute($dbID,$gid);
												}
											}
											else {
												$support->log_verbose("Not setting as a display_xref - Vega gene $gid($stripped_name) has a different name than the the MGI record ($xid)\n",1);
											next GENE;
											}
										}

										#hack to not set MGI as display xref unless it's a 'proper' name
										if ($support->param('xrefformat') eq 'mgi') {
											next GENE unless ($xid =~ /^[A-Z][a-z]{2}/);
										}

										#if there is a prefixed name then make a new DBEntry and set as display_xref
										if ($prefix) {
											my $new_xid = $prefix.':'.$stripped_name;
											my $new_dbentry = Bio::EnsEMBL::DBEntry->new(
																						 -primary_id => $pid,
																						 -display_id => $new_xid,
																						 -version    => 1,
																						 -release    => 1,
																						 -dbname     => $extdb,
																						 -info_text  => 'vega_source_prefix',
																						);
											my $new_dbID = $ea->store($new_dbentry, $gid, 'gene');
											$sth_display_xref->execute($new_dbID,$gid);
											$support->log("updated display_xref ($new_xid) using new $extdb xref ($new_dbID).\n",1);
										}

										#if the non-prefixed name matches the dbentry name then store it as display_xref
										elsif ($stripped_name eq $xid) {
											$sth_display_xref->execute($dbID,$gid);
											$support->log("updated display xref ($xid) using preexisting $extdb xref ($dbID).\n",1);
										}

										#warn if something has gone horribly wrong
										else {
											$support->log_warning("Expected name for display_xref ($xid) doesn't match Vega ($name), not updating Vega display_xref for gene $gsi ($gid).\n",1);
										}
									}	
								}
								else {
									$support->log_warning("No dbID for gene $gid, pid $pid, display_id $xid, dbname $extdb.\n", 1);
								}
							}						
							else {
								$support->log("Would store $extdb xref $xid (pid = $pid) for gene $gid.\n", 1);
								
							}
						}
					}
				}
			}
		# no match for some reason (log why)
		}
		if (!$xref_found) {
			foreach my $lc_name (@lc_names) {
				if ($lcmap->{$lc_name}) {
					# possible case error
					$support->log_warning("Possible case error for $gene_name -- ".
									  join(',',(@{ $lcmap->{$lc_name} }))."\n", 1);
					if ($support->param('mismatch')) {
						# fix case mismatch in db
						$support->log("Fixing case mismatch...\n", 1);
						my $new_name = $prefix ? $prefix.':'.$lcmap->{$lc_name}->[0] : $lcmap->{$lc_name}->[0];
#						$sth_case1->execute($new_name, $gene_name);
						$sth_case2->execute($new_name, $gene_name);
					}
					$warnings{'wrong_case'}++;
					next GENE;
				}
			}
			if ($gene_name =~ /^\w+\.\d+$/ || $gene_name =~ /^\w+\-\w+\.\d+$/) {
				# probably a clone-based genename - ok
				$support->log("No match for $gene_name (but has clonename based name).\n", 1);
				$warnings{'nomatch_clone'}++;
			}
			else {
				# other genes without a match
				$support->log_warning("No match for $gene_name (@syn_names).\n", 1);
				$warnings{'nomatch_missing'}++;
			}
		}
	}

    # log stats
    $support->log("\nProcessed $gnum genes (of ".scalar @$genes." on this chromosome).\n");
    $support->log("OK:\n");
    foreach my $extdb (sort keys %stats) {
        $support->log("$extdb $stats{$extdb}.\n", 1);
    }
    $support->log("WARNINGS:\n");
    $support->log("Genes with possible case mismatch: $warnings{wrong_case}.\n", 1);
    $support->log("Genes with apparently clonename based names: $warnings{nomatch_clone}.\n", 1);
    $support->log("Other genes without match: $warnings{nomatch_missing}.\n", 1);
    $support->log_stamped("Done with chromosome $chr.\n\n");

	$overall_stats{$chr} = \%stats;
}

#create a summary
my %report;
foreach my $chr_name (keys %overall_stats) {
	foreach my $extdb (sort keys %{$overall_stats{$chr_name}}) {
		$report{$extdb} += $overall_stats{$chr_name}->{$extdb};
	}
}

$support->log("\nSummary of xrefs assigned
-------------------------\n\n");
foreach my $extdb (keys %report) {
	$support->log("$extdb provides $report{$extdb} xrefs\n",1);
}


# finish log
$support->finish_log;


### end main ###


=head2 parse_ensdb

  Arg[1]      : Hashref $ens_xrefs - keys: gene names, values: hashref (extDB =>
                extID)
  Arg[2]      : Hashref $lcmap - keys: lowercase gene names, values: list of
                gene names (with case preserved)
  Arg[3]      : string - species short name
  Example     : &parse_ensdb($ens_xrefs, $lcmap);
                foreach my $gene (keys %$ens_xrefs) {
                    foreach my $extdb (keys %{ $ens_xrefs->{$gene} }) {
                        print "DB $extdb, extID ".$ens_xrefs->{$gene}->{$extdb}."\n";
                    }
                }
  Description : Parses xrefs from an E! core database
  Return type : none
  Exceptions  : thrown if database can't be read
  Caller      : internal

=cut

sub parse_ensdb {
	my ($ens_xrefs, $lcmap, $sp) = @_;
	$dba = $support->get_database('ensembl', 'ensembl');
	my $e_dbname = $support->param('ensembldbname');
    $support->log_stamped("Retrieving xrefs from $e_dbname...\n", 1);
	my $sa = $dba->get_SliceAdaptor();
	foreach my $chr ( @{$sa->fetch_all('chromosome')} ) {
		my $chr_name = $chr->seq_region_name;
		$support->log("Looking at chromosome $chr_name\n");
		foreach my $gene ( @{$chr->get_all_Genes} ) {
			my $stable_id = $gene->stable_id;
			my $disp_xref = $gene->display_xref;
			my $gene_name;
			if ($disp_xref) {
				$gene_name = $disp_xref->display_id;
			}
			next unless ($gene_name);
			no strict 'refs';
			if (exists($ens_xrefs->{$gene_name})) {
				$support->log_verbose("Ensembl gene $gene_name not unique, deleting stable id $stable_id\n",1);
				delete($ens_xrefs->{$gene_name}{'Ens_'.$sp.'_gene'});
			}
			else {
				$support->log_verbose("Storing Ensembl stable_id $stable_id\n",1);
				$ens_xrefs->{$gene_name}{'Ens_'.$sp.'_gene'} = $stable_id;
				#store lower case gene name for curation
				push @{ $lcmap->{lc($gene_name)} }, $gene_name;
			}
			my @db_entries = @{$gene->get_all_DBLinks()};
			foreach my $dbe (@db_entries){
				my $db_name = $dbe->dbname;
				my $dbid = $dbe->display_id(); 
				my $pid = $dbe->primary_id();
				push @{$ens_xrefs->{$gene_name}{$db_name}} , $dbid.'||'.$pid;
			}
			#parse display_xref to ensure it's added
			my $disp_dbname = $disp_xref->dbname;
			my $disp_dbid = $disp_xref->display_id();
			my $disp_pid  = $disp_xref->primary_id();
			$ens_xrefs->{$gene_name}{$disp_dbname} = $disp_dbid.'||'.$disp_pid;
			$support->log_verbose("Storing display_xref for gene $stable_id (db = $disp_dbname, display name = $disp_dbid, pid = $disp_pid\n");
		}
	}
}


=head2 parse_hugo

  Arg[1]      : Hashref $xrefs - keys: gene names, values: hashref (extDB =>
                extID)
  Arg[2]      : Hashref $lcmap - keys: lowercase gene names, values: list of
                gene names (with case preserved)
  Example     : &parse_hugo($xrefs, $lcmap);
                foreach my $gene (keys %$xrefs) {
                    foreach my $extdb (keys %{ $xrefs->{$gene} }) {
                        print "DB $extdb, extID ".$xrefs->{$gene}->{$extdb}."\n";
                    }
                }
  Description : Parses a nomeid file from HUGO.
  Return type : none
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_hugo {
    my ($xrefs, $lcmap) = @_;
	$support->log_stamped("Hugo...\n", 1);
    # read nomeid input file from HUGO
    open (NOM, '<', $support->param('hugofile')) or $support->throw(
        "Couldn't open ".$support->param('hugofile')." for reading: $!\n");
    # read header (containing external db names)
    my $line = <NOM>;
    chomp $line;
    my %convhash = (
        'HGNC ID'                      => 'HUGO',
        'UniProt ID (mapped data)'     => 'Uniprot/SWISSPROT',
        'RefSeq (mapped data)'         => 'RefSeq_dna',
        'Entrez Gene ID (mapped data)' => 'EntrezGene',
        'OMIM ID (mapped data)'        => 'MIM_GENE',
    );
    my @fieldnames = map { $convhash{$_} || $_ } split /\t/, $line;

    my $num_fields = scalar(@fieldnames);

    # parse input into datastructure
    my $alt_symbols;
    my %stats = (
        total           => 0,
        ok              => 0,
        missing_symbol  => 0,
    );
    while (<NOM>) {
        $stats{'total'}++;
        chomp;
        my @fields = split /\t/, $_, -1;

        # sanity checks
        if (scalar(@fields) != $num_fields) {
            $support->log("Wrong number of fields (got " . scalar(@fields) .
                          ", should be $num_fields).\n", 2);
            $support->log("Entry:\n$_\n", 2);
            $support->log("Aborting.\n", 2);
            die("Inconsistent field numbers in HUGO file. Aborting.\n");
        }
		my $hugo_pid = $fields[0];
        my $gene_name = $fields[1];

		#remove 'withdrawn' part of tyhe label since we still want to use these
		if ($gene_name =~ /\w+~withdrawn/) {
			$gene_name =~ s/(\w+)~withdrawn/$1/;
		}
		#parse input data
		my $i = 0;
        foreach my $fieldname (@fieldnames) {
			push @{$xrefs->{$gene_name}->{$fieldname}} , $fields[$i];
            $i++;
        }
		#modify HUGO entry for xref storing
		$xrefs->{$gene_name}->{'HUGO'}[0] = $gene_name .'||'. $hugo_pid;
		#modify EntrezGene Entry
		$xrefs->{$gene_name}->{'EntrezGene'}[0] .= '||'. $xrefs->{$gene_name}->{'EntrezGene'};
		#modify RefSeq_dna
		$xrefs->{$gene_name}->{'RefSeq_dna'}[0] .= '||'. $xrefs->{$gene_name}->{'RefSeq_dna'};
		#modify MIM Entry
		$xrefs->{$gene_name}->{'MIM_GENE'}[0] .=  '||'. $xrefs->{$gene_name}->{'MIM_GENE'};
		#modify Uniprot/SWISSPROT
		$xrefs->{$gene_name}->{'Uniprot/SWISSPROT'}[0] .= '||'.$xrefs->{$gene_name}->{'Uniprot/SWISSPROT'};

		#store lowercase name for matching
		push @{ $lcmap->{lc($gene_name)} }, $gene_name;

        # try to use previous symbols and aliases as a fallback as well
        my @previous = split(",", $fields[5]);
        my @aliases = split(",", $fields[7]);
        foreach my $alt_symbol (@previous, @aliases) {
           # trim whitespace
            $alt_symbol =~ s/^\s+//;
            $alt_symbol =~ s/\s+$//;
			
			#store details for alias
			unless ($alt_symbol eq $gene_name) {
				$xrefs->{$alt_symbol} = {%{$xrefs->{$gene_name}}};
				#make sure alias name and real pid are used
				$xrefs->{$alt_symbol}->{'HUGO'} = $alt_symbol .'||'. $hugo_pid;
				#add a remark for reporting
				$xrefs->{$alt_symbol}{'Is an alias for'} = $gene_name;
			}
        }
    }
    close(NOM);

    $support->log_stamped("Done processing ".$stats{'total'}." entries.\n\n", 1);
}

=head2 parse_mgivega

  Arg[1]      : Hashref $xrefs - keys: gene names, values: hashref (extDB =>
                extID) 
  Arg[2]      : Hashref $lcmap - keys: lowercase gene names, values: list of
                gene names (with case preserved)
  Example     : &parse_mgi($xrefs, $lcmap);
  Description : Parses a specific rtf file from MGI. Used to add MarkerSymbol, Swissprot
                RefSeq and EntrezGene xrefs to Vega genes
  Return type : none
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_mgivega {
   my ($xrefs) = @_;
   $support->log_stamped("MGI...\n", 1);

   # read input file
   my $mgivegafile = $support->param('mgivegafile');
   open(MGIV, "< $mgivegafile")
	   or $support->throw("Couldn't open $mgivegafile for reading: $!\n");

   #parse input file
   while (<MGIV>) {
	   chomp;
	   my @fields = split /\t/;
	   my $pid = $fields[0];
	   my $markersymbol = $fields[1];
       my $desc = $fields[2];
       my $vegaID = $fields[5];
	   if ( exists($xrefs->{$vegaID}) ) {
		   $support->log_warning("$vegaID found more than once in MGI file\n");
	   }
	   push @{$xrefs->{$vegaID}{'MarkerSymbol'}}, $markersymbol . '||' . $pid;
   }
}

=head2 parse_mgi

  Arg[1]      : Hashref $xrefs - keys: gene names, values: hashref (extDB =>
                extID) 
  Arg[2]      : Hashref $lcmap - keys: lowercase gene names, values: list of
                gene names (with case preserved)
  Example     : &parse_mgi($xrefs, $lcmap);
  Description : Parses a specific rtf file from MGI. Used to add MarkerSymbol, Swissprot
                RefSeq and EntrezGene xrefs to Vega genes
  Return type : none
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_mgi {
   my ($xrefs, $lcmap) = @_;
   $support->log_stamped("MGI...\n", 1);

   # read input file
   my $mgifile = $support->param('mgifile');
   open(MGI, "< $mgifile")
	   or $support->throw("Couldn't open $mgifile for reading: $!\n");

   #parse input file
   my %types;
   while (<MGI>) {
	   my @fields = split /\t/;

	   # MGI record contains all sorts of entries as well as genes, eg markers.
       # There doesn't seem to be any way of distinguishing between them, but could
	   # skip all apparent non-gene entries in the input file ?
#	   next if $fields[2] =~ /^RIKEN cDNA/;
#	   next if $fields[2] =~ /^DNA segment/;

	   #add mgi entries
	   my $gene_name = $fields[1];
	   my $mgi_pid = $fields[0];
	   $xrefs->{$gene_name}->{'MarkerSymbol'} = [ $gene_name .'||'. $mgi_pid ];

	   #add refseq dna entries
	   my $refseqs = $fields[16];
	   my (@ids) = split ',',$refseqs; 
	   foreach my $id (@ids) {
		   if ($id =~ /^NM_|NR_/) {
			   push @{$xrefs->{$gene_name}->{'RefSeq_dna'}}, $id.'||'.$id ;
		   }
		   if ($id =~ /^NP_/) {
			   push @{$xrefs->{$gene_name}->{'RefSeq_peptide'}}, $id.'||'.$id ;
		   }
	   }

	   #add swissprot entry
	   my $swissptrots = $fields[17];
	   my ($first_id) = split ',',$swissptrots; 
	   $xrefs->{$gene_name}->{'Uniprot/SWISSPROT'} = [ $first_id .'||'.$first_id ];

	   #add entrezgene entry
	   my $entrezgenes = $fields[20];
	   ($first_id) = split ',',$entrezgenes; 
	   $xrefs->{$gene_name}->{'EntrezGene'} = [ $first_id .'||'. $first_id ];

	   #add ensembl xrefs
	   my $ensid = $fields[9];
	   if ( $ensid =~ /^ENSMUSG/ ) {
		   $xrefs->{$gene_name}->{'Ens_Mm_gene'} = [ $ensid .'||'. $ensid ];
	   }
	   elsif ($ensid) {
		   $support->log_warning("Gene $gene_name from MGI has a non-mouse Ensembl ID ($ensid)\n");
	   }

	   #store lower case name to catch case mismatches
	   push @{ $lcmap->{lc($gene_name)} }, $gene_name;
   }

}



=head2 parse_tcag

=cut

sub parse_tcag {
    my ($xrefs, $lcmap) = @_;
	$support->log_stamped("TCAG...\n", 1);

    # read input file
    my $tcagfile = $support->param('tcagfile');
    my $fh_expr;
    if($tcagfile =~ /\.gz$/) {
        $fh_expr = "gzip -d -c $tcagfile |";
    } else {
        $fh_expr = "< $tcagfile";
    }
    open(TCAG, $fh_expr)
        or $support->throw("Couldn't open $tcagfile for reading: $!\n");
	
	#parse input file
	while (<TCAG>) {
		my @fields = split /\t/;
		next unless  $fields[8] =~ /Gene_ID/;
		my @details = split /;/, $fields[8];
 		my ($symbol) = $details[0] =~ /"(.+)"/;
		next if ($symbol =~ /transcript_variant/);
		my ($id)     = $details[2] =~ /"(.+)"/;

		#for debugging
		unless ($symbol && $id) {
			print "record with no symbol = $_\n";
			foreach my $f (@fields) {
				print " field = $f\n";
			}
			foreach my $d (@details) {
				print "  det = $d\n";
			}
		}

		push @{$xrefs->{$symbol}->{'TCAG'}} , $id.'||'.$id;
		push @{ $lcmap->{lc($symbol)} }, $symbol;
	}
}

=head2 parse_imgt

=cut

sub parse_imgt {
    my ($xrefs, $lcmap) = @_;
	$support->log_stamped("IMGT...\n", 1);
    # read input file from IMGT
    open (IMGT, '<', $support->param('imgtfile')) or $support->throw(
        "Couldn't open ".$support->param('imgtfile')." for reading: $!\n");
    # read header
    my $line = <IMGT>;
    chomp $line;
	my @fieldnames = split /\t/, $line;
	while (<IMGT>) {
		my @fields = split /\t/, $_;
		#skip multiple header lines
		next if ($fields[0] eq $fieldnames[1]);
		my $xid  = $fields[0];
		my $gsi = $fields[2];
		my $pid = $xid;
		$pid =~ s/\*/_/; 
		push @{$xrefs->{$gsi}->{'IMGT'}} , $xid.'||'.$pid;
	}
}


=head2 parse_locuslink

  Arg[1]      : Hashref $xrefs - keys: gene names, values: hashref (extDB =>
                extID)
  Arg[2]      : Hashref $lcmap - keys: lowercase gene names, values: list of
                gene names (with case preserved)
  Example     : &parse_locuslink($xrefs, $lcmap);
                foreach my $gene (keys %$xrefs) {
                    foreach my $extdb (keys %{ $xrefs->{$gene} }) {
                        print "DB $extdb, extID ".$xrefs->{$gene}->{$extdb}."\n";
                    }
                }
  Description : Parses a LL_tmpl file from LocusLink. File can optionally be
                gzipped.
  Return type : none
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_locuslink {
    my ($xrefs, $lcmap) = @_;

    $support->log_stamped("LocusLink...\n", 1);

    # read LL_tmpl input file
    my $xreffile = $support->param('locuslinkfile');
    my $fh_expr;
    if($xreffile =~ /\.gz$/) {
        $fh_expr = "gzip -d -c $xreffile |";
    } else {
        $fh_expr = "< $xreffile";
    }
    open(LL, $fh_expr)
        or $support->throw("Couldn't open $xreffile for reading: $!\n");

    # parse input into datastructure
    my %stats = (
        total           => 1,
        ok              => 0,
        wrong_org       => 0,
        missing_org     => 0,
        missing_symbol  => 0,
    );
    my $locus_id;
    my $flag_found = 1;
    my $flag_org = 0;
    my $nm = '';
    my $mgi = '';
    my $gene_name = '';

    while(<LL>){
        if (/^NM: (.*)/) {
            $nm = $1;
            $nm =~ s/\|.*//;
        } elsif (/http:\/\/www.informatics.jax.org\/searches\/accession_report.cgi\?id=MGI:(\d+)/) {
            $mgi = $1;
        } elsif (/^ORGANISM: (.*)/) {
            my $org = $1;
            if ($species eq $org) {
                $flag_org = 1;
            } else {
                $flag_org = 2;
            }
        } elsif ((/^OFFICIAL_SYMBOL: (\w+)/) || (/^PREFERRED_SYMBOL: (\w+)/)) {
            if ($flag_org == 1) {
                $gene_name = $1;
                push @{ $lcmap->{lc($gene_name)} }, $gene_name;
                $stats{'ok'}++;
            } elsif($flag_org == 2) {
                # wrong organism
                $stats{'wrong_org'}++;
            } elsif($flag_org == 0) {
                # missing organism
                $support->log_warning("ORGANISM not found for $locus_id.\n", 2);
                $stats{'missing_org'}++;
            }
            $flag_found = 1;
        } elsif (/\>\>(\d+)/) {
            if ($gene_name) {
                $xrefs->{$gene_name} = {
                    RefSeq_dna      => $nm,
                    EntrezGene      => $locus_id,
                    MarkerSymbol    => $mgi,
                };
                $support->log_verbose("Gene $gene_name: EntrezGene ID $locus_id\n", 2);
            }
            unless ($flag_found) {
                $stats{'missing_symbol'}++;
                #$support->log_warning("OFFICIAL_SYMBOL for $locus_id not found.\n", 1);
            }
            $locus_id = $1;
            $stats{'total'}++;

            $nm = '';
            $mgi = '';
            $gene_name = '';
            $flag_found = 0;
            $flag_org = 0;
        }
    }
    close(LL);

    # log stats
    $support->log_stamped("Done processing ".$stats{'total'}." entries.\n", 1);
    $support->log("OK: ".$stats{'ok'}.".\n", 1);
    $support->log("SKIPPED:\n", 1);
    $support->log("Not \'$species\': ".$stats{'wrong_org'}.".\n", 2);
    $support->log("No organism label: ".$stats{'missing_org'}.".\n", 2);
    $support->log("No official symbol: ".$stats{'missing_symbol'}.".\n\n", 2);
}

=head2 parse_refseq

  Arg[1]      : Hashref $xrefs - keys: gene names, values: hashref (extDB =>
                extID)
  Arg[2]      : Hashref $lcmap - keys: lowercase gene names, values: list of
                gene names (with case preserved)
  Example     : &parse_refseq($xrefs, $lcmap);
                foreach my $gene (keys %$xrefs) {
                    foreach my $extdb (keys %{ $xrefs->{$gene} }) {
                        print "DB $extdb, extID ".$xrefs->{$gene}->{$extdb}."\n";
                    }
                }
  Description : Parses a xref file from RefSeq in Genbank format
  Return type : none
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_refseq {
    my ($xrefs, $lcmap) = @_;

    $support->log_stamped("Refseq...\n", 1);
    
    # read input from refseq file (genbank format)
    my $in = Bio::SeqIO->new(
        -file => $support->param('refseqfile'),
        -format => 'genbank'
    );

    # parse input into datastructure
    my %stats = (
        total           => 0,
        ok              => 0,
        not_nm          => 0,
        missing_symbol  => 0,
    );
    my $found = 0;

    while (my $seq = $in->next_seq) {
        my $id = $seq->id;
        $stats{'total'}++;

        # only use NM identifiers
        unless ($id =~ /NM_.*/) {
            $stats{'not_nm'}++;
            next;
        }

        foreach my $f ($seq->get_SeqFeatures) {
            if ($f->has_tag('gene')) {
                foreach my $gene_name ($f->get_tag_values('gene')) {
                    $xrefs->{$gene_name} = { RefSeq_dna => $id };
                    push @{ $lcmap->{lc($gene_name)} }, $gene_name;
                    $stats{'genenames'}++;
                    $found = 1;
                }
            }
        }

        $stats{'missing_symbol'}++ unless $found;
        $found = 0;
    }

    # exit if no entries were found
    unless ($stats{'total'}) {
        $support->throw("No entries found (have you used the right input format?). Aborting...\n");
    }

    # log stats
    $support->log_stamped("Done processing ".$stats{'total'}." entries.\n", 1);
    $support->log("OK:\n", 1);
    $support->log("Found ".$stats{'genenames'}." gene names (".scalar(keys(%$xrefs))." unique).\n", 2);
    $support->log("SKIPPED:\n", 1);
    $support->log("No NM identifier: ".$stats{'not_nm'}.".\n", 2);
    $support->log("No gene symbol: ".$stats{'missing_symbol'}.".\n\n", 2);
}

