#!/usr/local/bin/perl

=head1 NAME

add_external_xrefs.pl - adds xrefs to external databases from various types
of input files

=head1 SYNOPSIS

    add_external_xrefs.pl [options]

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

    Specific options:
        --chromosomes, --chr=LIST           only process LIST chromosomes
        --gene_stable_id, --gsi=LIST|FILE   only process LIST gene_stable_ids
                                            (or read list from FILE)
        --xreffile FILE                     read input from FILE
        --xrefformat FORMAT                 input file format FORMAT
                                            (hugo|locuslink|refseq)

=head1 DESCRIPTION

This script parses input files from various sources (HUGO, LocusLink, RefSeq)
and adds xrefs to the databases covered by the respective input source. If
appropriate, display name of genes is set accordingly.

Currently, these input formats are supported:
    
    hugo        => http://www.gene.ucl.ac.uk/nomenclature/nomeids.txt
    locuslink   => ftp://ftp.ncbi.nih.gov/refseq/LocusLink/LL_tmpl.gz
    refseq      => ftp://ftp.ncbi.nih.gov/genomes/__SPECIES__/RNA/rna.gkb

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

my $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
    'chromosomes|chr=s@',
    'gene_stable_id|gsi=s@',
    'xreffile=s',
    'xrefformat=s',
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

$support->comma_to_list('chromosomes');
$support->list_or_file('gene_stable_id');

# ask user to confirm parameters to proceed
$support->confirm_params;

# make sure add_vega_xrefs.pl has been run
print "This script must run after add_vega_xrefs.pl. Have you run it?\n";
$support->user_confirm;

# get log filehandle and print heading and parameters to logfile
$support->log_filehandle('>>');
$support->log($support->init_log);

# connect to database and get adaptors (caching features on one slice only)
# get an ensembl database for better performance (no otter tables are needed)
my $dba = $support->get_database('ensembl');
$Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::SLICE_FEATURE_CACHE_SIZE = 1;
my $sa = $dba->get_SliceAdaptor();
my $ga = $dba->get_GeneAdaptor();
my $ea = $dba->get_DBEntryAdaptor();
# statement handle for display_xref_id update
my $sth = $dba->dbc->prepare("UPDATE gene SET display_xref_id=? WHERE gene_id=?");

my @gene_stable_ids = $support->param('gene_stable_id');
my %gene_stable_ids = map { $_, 1 } @gene_stable_ids;
my $chr_length = $support->get_chrlength($dba);
my @chr_sorted = $support->sort_chromosomes($chr_length);

# get species name
my $species = $support->get_species_scientific_name($dba);

# sanity checks
$support->check_required_params('xreffile', 'xrefformat');
my %allowed_formats = map { $_ => 1 } qw(hugo refseq locuslink);
unless ($allowed_formats{$support->param('xrefformat')}) {
    $support->throw("Invalid xrefformat ".$support->param('xrefformat').".\nAllowed formats: ".join(" ", keys %allowed_formats));
}

# parse input file
$support->log("Reading xref input file... ".$support->date_and_mem."\n");
my $parser = 'parse_'.$support->param('xrefformat');
no strict 'refs';
my ($xrefs, $lcmap) = &$parser($support);

# define what to do with each type of xref
my %primary = ( $support->param('xrefformat') => 1 );
my %extdb_def = (
    HUGO                => [ 'KNOWN', $primary{'hugo'} ],
    LocusLink           => [ 'KNOWNXREF', $primary{'locuslink'} ],
    RefSeq_dna          => [ 'KNOWN', 0 ],
    MIM                 => [ 'KNOWNXREF', 0 ],
    'Uniport/SWISSPROT' => [ 'KNOWN', 0 ],
);

# loop over chromosomes
$support->log("Looping over chromosomes: @chr_sorted\n\n");
foreach my $chr (@chr_sorted) {
    $support->log("> Chromosome $chr (".$chr_length->{$chr}
               ."bp). ".$support->date_and_mem."\n\n");
    
    # fetch genes from db
    $support->log("Fetching genes...\n");
    my $slice = $sa->fetch_by_region('chromosome', $chr);
    my $genes = $ga->fetch_all_by_Slice($slice);
    $support->log("Done fetching ".scalar @$genes." genes. " .
                   $support->date_and_mem."\n\n");

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
        my $gene_name = $gene->display_xref->display_id;
        my $lc_gene_name = lc($gene_name);

        # filter to user-specified gene_stable_ids
        if (scalar(@gene_stable_ids)){
            next unless $gene_stable_ids{$gsi};
        }

        $gnum++;
        $support->log("Gene $gene_name ($gid, $gsi)...\n");

        # we have a match
        if ($xrefs->{$gene_name}) {
            foreach my $extdb (keys %extdb_def) {
                if (my $xid = $xrefs->{$gene_name}->{$extdb}) {
                    $stats{$extdb}++;
                    my $display_id;
                    if ($extdb_def{$extdb}->[1]) {
                        $display_id = $gene_name;
                    } else {
                        $display_id = $xid;
                    }
                    my $dbentry = Bio::EnsEMBL::DBEntry->new(
                            -primary_id => $xid,
                            -display_id => $display_id,
                            -version    => 1,
                            -release    => 1,
                            -dbname     => $extdb,
                    );
                    $dbentry->status($extdb_def{$extdb}->[0]);
                    $gene->add_DBEntry($dbentry);
                    if ($support->param('dry_run')) {
                        $support->log("Would store $extdb xref $xid for gene $gid.\n", 1);
                    } else {
                        $ea->store($dbentry, $gid, 'Gene');
                        if ($extdb_def{$extdb}->[1]) {
                            $sth->execute($dbentry->dbID, $gid);
                        }
                        $support->log("Stored $extdb xref $xid (dbID ".$dbentry->dbID.") for gene $gid.\n", 1);
                    }
                }
            }

        # no match for some reason (log why)
        } elsif ($lcmap->{$lc_gene_name}) {
            # possible case error
            $support->log_warning("Possible case error for $gene_name: ".
                join(',',(@{ $lcmap->{$lc_gene_name} })).".\n", 1);
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
    $support->log("Done with chromosome $chr. ".$support->date_and_mem."\n\n");
}

# finish log
$support->log($support->finish_log);


### end main ###


=head2 parse_hugo

  Arg[1]      : Bio::EnsEMBL::Utils::ConversionSupport object $support
  Example     : my ($xrefs, $lcmap) = &parse_hugo($support);
                foreach my $gene (keys %$xrefs) {
                    foreach my $extdb (keys %{ $xrefs->{$gene} }) {
                        print "DB $extdb, extID ".$xrefs->{$gene}->{$extdb}."\n";
                    }
                }
  Description : Parses a nomeid file from HUGO.
  Return type : list of hashrefs:
                1.  keys: gene_names
                    values: hashref (extDB => extID)
                2.  keys: lowercase gene_names
                    values: list of gene_names (with case preserved)
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_hugo {
    my $support = shift;

    # read nomeid input file from HUGO
    open (NOM, '<', $support->param('xreffile')) or $support->throw(
        "Couldn't open ".$support->param('xreffile')." for reading: $!\n");

    # read header (containing external db names)
    my $line = <NOM>;
    chomp $line;
    my %convhash = (
        'HGNC'          => 'HUGO',
        'SWISSPROT'     => 'Uniprot/SWISSPROT',
        'Ref Seq'       => 'RefSeq_dna',
        'Locus Link'    => 'LocusLink',
    );
    my @fieldnames = map { $convhash{$_} || $_ } split /\t/, $line;
    my $num_fields = scalar(@fieldnames);

    # parse input into datastructure
    my (%xrefs, %lcmap);
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
                          ", should be $num_fields).\n");
            $support->log("Entry:\n$_\n");
            $support->log("Aborting.\n");
            die("Inconsistent field numbers in HUGO file. Aborting.\n");
        }
        if (defined($xrefs{$fields[1]})) {
            $support->log("Duplicate lines for " . $fields[1] . ". Aborting.\n");
            die("Duplicate lines in HUGO file. Aborting.\n");
        }
    
        my $gene_name = $fields[1];
        @{ $xrefs{$gene_name} }{@fieldnames} = @fields;
        push @{ $lcmap{lc($gene_name)} }, $gene_name;
    }
    close(NOM);

    $support->log("Done processing ".$stats{'total'}." entries. ".$support->date_and_mem."\n\n");

    return (\%xrefs, \%lcmap);
}

=head2 parse_locuslink

  Arg[1]      : Bio::EnsEMBL::Utils::ConversionSupport object $support
  Example     : my ($xrefs, $lcmap) = &parse_locuslink($support);
                foreach my $gene (keys %$xrefs) {
                    foreach my $extdb (keys %{ $xrefs->{$gene} }) {
                        print "DB $extdb, extID ".$xrefs->{$gene}->{$extdb}."\n";
                    }
                }
  Description : Parses a LL_tmpl file from LocusLink. File can optionally be
                gzipped.
  Return type : list of hashrefs:
                1.  keys: gene_names
                    values: hashref (extDB => extID)
                2.  keys: lowercase gene_names
                    values: list of gene_names (with case preserved)
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_locuslink {
    my $support = shift;

    # read LL_tmpl input file
    my $xreffile = $support->param('xreffile');
    my $fh_expr;
    if($xreffile =~ /\.gz$/) {
        $fh_expr = "gzip -d -c $xreffile |";
    } else {
        $fh_expr = "< $xreffile";
    }
    open(LL, $fh_expr)
        or $support->throw("Couldn't open $xreffile for reading: $!\n");

    # parse input into datastructure
    my (%xrefs, %lcmap);
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
    my $gene_name = '';

    while(<LL>){
        if (/^NM: (.*)/) {
            $nm = $1;
            $nm =~ s/\|.*//;
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
                push @{ $lcmap{lc($gene_name)} }, $gene_name;
                $stats{'ok'}++;
            } elsif($flag_org == 2) {
                # wrong organism
                $stats{'wrong_org'}++;
            } elsif($flag_org == 0) {
                # missing organism
                $support->log_warning("ORGANISM not found for $locus_id.\n", 1);
                $stats{'missing_org'}++;
            }
            $flag_found = 1;
        } elsif (/\>\>(\d+)/) {
            if ($gene_name) {
                $xrefs{$gene_name} = {
                    RefSeq_dna  => $nm,
                    LocusLink   => $locus_id,
                };
                $support->log("Gene $gene_name: LocusLink ID $locus_id\n", 1);
            }
            unless ($flag_found) {
                $stats{'missing_symbol'}++;
                #$support->log_warning("OFFICIAL_SYMBOL for $locus_id not found.\n", 1);
            }
            $locus_id = $1;
            $stats{'total'}++;

            $nm = '';
            $gene_name = '';
            $flag_found = 0;
            $flag_org = 0;
        }
    }
    close(LL);

    # log stats
    $support->log("Done processing ".$stats{'total'}." entries. ".$support->date_and_mem."\n");
    $support->log("OK: ".$stats{'ok'}.".\n");
    $support->log("SKIPPED:\n");
    $support->log("Not \'$species\': ".$stats{'wrong_org'}.".\n", 1);
    $support->log("No organism label: ".$stats{'missing_org'}.".\n", 1);
    $support->log("No official symbol: ".$stats{'missing_symbol'}.".\n\n", 1);

    return (\%xrefs, \%lcmap);
}

=head2 parse_refseq

  Arg[1]      : Bio::EnsEMBL::Utils::ConversionSupport object $support
  Example     : my ($xrefs, $lcmap) = &parse_refseq($support);
                foreach my $gene (keys %$xrefs) {
                    foreach my $extdb (keys %{ $xrefs->{$gene} }) {
                        print "DB $extdb, extID ".$xrefs->{$gene}->{$extdb}."\n";
                    }
                }
  Description : Parses a xref file from RefSeq in Genbank format
  Return type : list of hashrefs:
                1.  keys: gene_names
                    values: hashref (extDB => extID)
                2.  keys: lowercase gene_names
                    values: list of gene_names (with case preserved)
  Exceptions  : thrown if input file can't be read
  Caller      : internal

=cut

sub parse_refseq {
    my $support = shift;

    # read input from refseq file (genbank format)
    my $in = Bio::SeqIO->new(
        -file => $support->param('xreffile'),
        -format => 'genbank'
    );

    # parse input into datastructure
    my (%xrefs, %lcmap);
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
                    $xrefs{$gene_name} = { RefSeq_dna => $id };
                    push @{ $lcmap{lc($gene_name)} }, $gene_name;
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
    $support->log("Done processing ".$stats{'total'}." entries. ".$support->date_and_mem."\n");
    $support->log("OK:\n");
    $support->log("Found ".$stats{'genenames'}." gene names (".scalar(keys(%xrefs))." unique).\n", 1);
    $support->log("SKIPPED:\n");
    $support->log("No NM identifier: ".$stats{'not_nm'}.".\n", 1);
    $support->log("No gene symbol: ".$stats{'missing_symbol'}.".\n\n", 1);

    return (\%xrefs, \%lcmap);
}

