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
                                        (hugo|locuslink|refseq)
    --hugofile=FILE                     read Hugo input from FILE
    --locuslinkfile=FILE                read LocusLink input from FILE
    --refseqfile=FILE                   read Refseq input from FILE
    --mismatch                          correct case mismatches in the db
                                        (NOTE: this option overrides
                                            dry_run!)

=head1 DESCRIPTION

This script parses input files from various sources (HUGO, LocusLink, RefSeq)
and adds xrefs to the databases covered by the respective input source. If
appropriate, display name of genes is set accordingly.

Currently, these input formats are supported:
    
    hugo        => http://www.gene.ucl.ac.uk/nomenclature/data/gdlw_index.html
                   ('All data' in text format)
    locuslink   => ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz
    refseq      => ftp://ftp.ncbi.nih.gov/genomes/__SPECIES__/RNA/rna.gkb

For human, a combination of locuslink and hugo can be used by specifying both
as source (NOTE: order matters, locuslink has to be read first!)

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>
Original code by Tim Hubbard <th@sanger.ac.uk>

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
use Bio::SeqIO::genbank;

$| = 1;

our $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'chromosomes|chr=s@',
    'gene_stable_id|gsi=s@',
    'xrefformat=s@',
    'hugofile=s',
    'locuslinkfile=s',
    'refseqfile=s',
    'mismatch',
);
$support->allowed_params(
    $support->get_common_params,
    'chromosomes',
    'gene_stable_id',
    'xrefformat',
    'hugofile',
    'locuslinkfile',
    'refseqfile',
    'mismatch',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');
$support->comma_to_list('xrefformat');
$support->list_or_file('gene_stable_id');

# ask user to confirm parameters to proceed
$support->confirm_params;

# make sure add_vega_xrefs.pl has been run
exit unless $support->user_proceed("This script must run after add_vega_xrefs.pl. Have you run it?");

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
$Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::SLICE_FEATURE_CACHE_SIZE = 1;
my $sa = $dba->get_SliceAdaptor();
my $ga = $dba->get_GeneAdaptor();
my $ea = $dba->get_DBEntryAdaptor();

# statement handle for display_xref_id update
my $sth_display_xref = $dba->dbc->prepare("UPDATE gene SET display_xref_id=? WHERE gene_id=?");

# statement handles for fixing case errors
my $sth_case1 = $dba->dbc->prepare("UPDATE gene_name set name = ? WHERE name = ?");
my $sth_case2 = $dba->dbc->prepare("UPDATE xref set display_label = ? WHERE display_label = ?");

my @gene_stable_ids = $support->param('gene_stable_id');
my %gene_stable_ids = map { $_, 1 } @gene_stable_ids;
my $chr_length = $support->get_chrlength($dba);
my @chr_sorted = $support->sort_chromosomes($chr_length);

# get species name
my $species = $support->get_species_scientific_name($dba);

# sanity checks
$support->check_required_params('xrefformat');
my %allowed_formats = map { $_ => 1 } qw(hugo refseq locuslink);

# parse input file
$support->log_stamped("Reading xref input files...\n");
no strict 'refs';
my %primary;
my $xrefs = {};
my $lcmap = {};
foreach my $format ($support->param('xrefformat')) {
    my $parser = "parse_$format";
    &$parser($xrefs, $lcmap);

    # set as primary xref (will be used as display_xref)
    $primary{$format} = 1;
}
use strict 'refs';
$support->log_stamped("Done.\n\n");

# define what to do with each type of xref
my %extdb_def = (
    HUGO                => [ 'KNOWN', $primary{'hugo'} ],
    EntrezGene          => [ 'KNOWNXREF', $primary{'hugo'} || $primary{'locuslink'} ],
    MarkerSymbol        => [ 'KNOWNXREF', $primary{'locuslink'} ],
    RefSeq_dna          => [ 'KNOWN', 0 ],
    MIM                 => [ 'KNOWNXREF', 0 ],
    'Uniprot/SWISSPROT' => [ 'KNOWN', 0 ],
);

# loop over chromosomes
$support->log("Looping over chromosomes: @chr_sorted\n\n");
my $seen_xrefs;
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
    foreach my $gene (@$genes) {
        my $gsi = $gene->stable_id;
        my $gid = $gene->dbID;
        
        # catch missing display_xrefs here!!
        my $disp_xref = $gene->display_xref;
        my $gene_name;
        if ($disp_xref) {
            $gene_name = $disp_xref->display_id;
        } else {
            $support->log_warning("No display_xref found for gene $gid ($gsi). Skipping.\n");
            next;
        }

        # strip prefixes from gene names
        (my $stripped_name = $gene_name) =~ s/.*:(.*)/$1/;

        my $lc_gene_name = lc($gene_name);

        # filter to user-specified gene_stable_ids
        if (scalar(@gene_stable_ids)){
            next unless $gene_stable_ids{$gsi};
        }

        $gnum++;
        $support->log("Gene $gene_name ($gid, $gsi)...\n");

        # we have a match
        if ($xrefs->{$gene_name}) {
            # the sort is important so that MarkerSymbol superseeds EntrezGene
            # when setting display_xrefs
            foreach my $extdb (sort keys %extdb_def) {
                if (my $xid = $xrefs->{$gene_name}->{$extdb}) {
                    $stats{$extdb}++;
                    my $display_id;
                    if ($extdb_def{$extdb}->[1]) {
                        $display_id = $gene_name;
                    } else {
                        $display_id = $xid;
                    }

                    # check if we already had an xref with the sme primary_id
                    # and dbname, but different display_id; if so, bump up the
                    # version to force the xref to be stored
                    $seen_xrefs->{"$extdb:$xid"}->{$display_id} = 1;
                    my $version = scalar(keys(%{ $seen_xrefs->{"$extdb:$xid"} }));

                    my $dbentry = Bio::EnsEMBL::DBEntry->new(
                            -primary_id => $xid,
                            -display_id => $display_id,
                            -version    => $version,
                            -release    => 1,
                            -dbname     => $extdb,
                    );
                    $dbentry->status($extdb_def{$extdb}->[0]);
                    $gene->add_DBEntry($dbentry);
                    if ($support->param('dry_run')) {
                        $support->log("Would store $extdb xref $xid for gene $gid.\n", 1);
                    } else {
                        $ea->store($dbentry, $gid, 'Gene');
                        my $dbID;
                        while ($dbID == 0) {
                            $dbID = $dbentry->dbID;
                        }
                        $support->log("Storing $extdb xref $xid (dbID $dbID) for gene $gid.\n", 1);
                        if ($extdb_def{$extdb}->[1]) {
                            $sth_display_xref->execute($dbID, $gid);
                            $support->log("Setting display_xref_id to $dbID.\n", 1);
                        }
                    }
                }
            }

        # we have a match with the stripped gene name
        } elsif ($xrefs->{$stripped_name}) {
            foreach my $extdb (sort keys %extdb_def) {
                if (my $xid = $xrefs->{$stripped_name}->{$extdb}) {
                    $stats{$extdb}++;
                    my $display_id;
                    if ($extdb_def{$extdb}->[1]) {
                        $display_id = $stripped_name;
                    } else {
                        $display_id = $xid;
                    }

                    $seen_xrefs->{"$extdb:$xid"}->{$display_id} = 1;
                    my $version = scalar(keys(%{ $seen_xrefs->{"$extdb:$xid"} }));

                    my $dbentry = Bio::EnsEMBL::DBEntry->new(
                            -primary_id => $xid,
                            -display_id => $display_id,
                            -version    => $version,
                            -release    => 1,
                            -dbname     => $extdb,
                    );
                    $dbentry->status($extdb_def{$extdb}->[0]);
                    $gene->add_DBEntry($dbentry);
                    if ($support->param('dry_run')) {
                        $support->log("Would store $extdb xref $xid for gene $gid.\n", 1);
                    } else {
                        $ea->store($dbentry, $gid, 'Gene');
                        my $dbID;
                        while ($dbID == 0) {
                            $dbID = $dbentry->dbID;
                        }
                        $support->log("Storing $extdb xref $xid (dbID $dbID) for gene $gid.\n", 1);
                        if ($extdb_def{$extdb}->[1]) {
                            $sth_display_xref->execute($dbID, $gid);
                            $support->log("Setting display_xref_id to $dbID.\n");
                        }
                    }
                }
            }

        # no match for some reason (log why)
        } elsif ($lcmap->{$lc_gene_name}) {
            # possible case error
            $support->log_warning("Possible case error for $gene_name: ".
                join(',',(@{ $lcmap->{$lc_gene_name} })).".\n", 1);
            if ($support->param('mismatch')) {
                # fix case mismatch in db
                $support->log("Fixing case mismatch...\n", 1);
                $sth_case1->execute($lcmap->{$lc_gene_name}->[0], $gene_name);
                $sth_case2->execute($lcmap->{$lc_gene_name}->[0], $gene_name);
            }
            $warnings{'wrong_case'}++;
        } elsif ($gene_name =~ /^\w+\.\d+$/ || $gene_name =~ /^\w+\-\w+\.\d+$/) {
            # probably a clone-based genename - ok
            $support->log("No match for $gene_name (but has clonename based name).\n", 1);
            $warnings{'nomatch_clone'}++;
        } else {
            # other genes without a match
            $support->log_warning("No match for $gene_name.\n", 1);
            $warnings{'nomatch_missing'}++;
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
}

# finish log
$support->finish_log;


### end main ###


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
        'HGNC ID'                   => 'HUGO',
        'UniProt ID (mapped data)'  => 'Uniprot/SWISSPROT',
        'RefSeq (mapped data)'      => 'RefSeq_dna',
        'Entrez Gene ID'            => 'EntrezGene',
        'OMIM ID (mapped data)'     => 'MIM',
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
    
        my $gene_name = $fields[1];

        # complement xrefs from locuslink run
        my $i = 0;
        foreach my $fieldname (@fieldnames) {
            $xrefs->{$gene_name}->{$fieldname} ||= $fields[$i];
            $i++;
        }
        
        push @{ $lcmap->{lc($gene_name)} }, $gene_name;

        # try to use previous symbols and aliases as a fallback as well
        my @previous = split(",", $fields[5]);
        my @aliases = split(",", $fields[7]);
        foreach my $alt_symbol (@previous, @aliases) {
            # trim whitespace
            $alt_symbol =~ s/^\s+//;
            $alt_symbol =~ s/\s+$//;

            # store alternative name mapping in temporary hash
            @{ $alt_symbols->{$alt_symbol} }{@fieldnames} = @fields;
            push @{ $lcmap->{lc($alt_symbol)} }, $alt_symbol;
        }
    }
    close(NOM);

    # now loop over alternative symbols and add them to %xref if the symbol
    # name doesn't exist yet
    foreach my $symbol (keys %$alt_symbols) {
        foreach my $fieldname (@fieldnames) {
            $xrefs->{$symbol}->{$fieldname} ||= $alt_symbols->{$symbol}->{$fieldname};
        }
    }

    $support->log_stamped("Done processing ".$stats{'total'}." entries.\n\n", 1);
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
        } elsif (/^OFFICIAL_SYMBOL: (\w+)/) {
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

