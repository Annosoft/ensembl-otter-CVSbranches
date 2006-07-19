#!/usr/local/bin/perl

=head1 NAME

Insert description into the analysis_description table

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
                                        (hugo|locuslink|refseq|ensemblxref)
    --hugofile=FILE                     read Hugo input from FILE
    --locuslinkfile=FILE                read LocusLink input from FILE
    --refseqfile=FILE                   read Refseq input from FILE
    --mismatch                          correct case mismatches in the db
                                        (NOTE: this option overrides
                                            dry_run!)
    --prune                             reset to the state before running this
                                        script (i.e. after running
                                        add_vega_xrefs.pl)

=head1 DESCRIPTION
This script inserts descriptions into the analysis_descriptions table, based on the logic name in the analysis table. The description text is hard-coded
in this script.

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Stephen Rice (sr7@sanger.ac.uk)

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
    
    #$SERVERROOT = "/ensweb/sr7";
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

$support->allowed_params(
    $support->get_common_params
    
);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}


# ask user to confirm parameters to proceed
$support->confirm_params;


#ask user if they have run the prerequisite script that tags CORF genes

if(! $support->user_proceed("IMPORTANT: proceed only if you have already run 'change_analysis_by_stable_id'.\nContinue?")){

  exit(1);

};



# get log filehandle and print heading and parameters to logfile
$support->init_log;


my %desc_mapper= (
'Fgenesh'  => 'This transcript was predicted by the Vega pipeline analysis system using Fgenesh, a HMM-based gene-finding program with an algorithm similar to that in Genscan (Burge, C. and Karlin, S. 1997. J. Mol. Biol. 268, 78-94) and Genie (Kulp et al. 1996. A generalized hidden Markov model for the recognition of human genes in DNA. Intell. Systems Mol. Biol. 4: 134-142; Salamov A.A. and Solovyev V.V. 2000. Genome Research 2000 Apr; 10(4): 516-522)',
'Genscan' => 'This transcript was predicted by the Vega pipeline analysis system using Genscan (Burge, C. and Karlin, S. 1997. Prediction of complete gene structures in human genomic DNA. J. Mol. Biol. 268, 78-94; Burge, C. 1998. Modeling dependencies in pre-mRNA splicing signals, <i>in</i> Salzberg, S., Searls, D. and Kasif, S., eds. Computational Methods in Molecular Biology, Elsevier Science, Amsterdam, 127-163)',
'otter' => 'Finished genomic sequence is analysed on a clone by clone basis using a combination of similarity searches against DNA and protein databases as well as a series of ab initio gene predictions (GENSCAN, Fgenes). In addition, comparative analysis using vertebrate datasets is used to aid novel gene discovery. The data gathered in these steps is then used to manually annotate the clone adding gene structures, descriptions and poly-A features. The annotation is based on supporting evidence only.',
'otter_corf' => 'In contrast to the annotation on the fully annotated chromosomes, which attempts to identify every transcript at a particular loci, annotation for the CORF project at the WTSI is less complete, often only concentrating on transcripts with CCDS identifiers. See curation method of Vega genes for a full description of the annotation process.',
'otter_igsf' => 'In contrast to the annotation on the fully annotated chromosomes, which attempts to identify every transcript at a particular loci, annotation for the Atlas project at the WTSI is less complete, often only annotating the whole transcripts with the longest CDS. See curation method of Vega genes for a full description of the annotation process.'
);


# connect to database and get adaptors (caching features on one slice only)
# get an ensembl database for better performance (no otter tables are needed)
my $dba = $support->get_database('ensembl');

my $dbh= $dba->dbc->db_handle;


$support->log("Deleting the current contents of the 'analysis_description' table.\n",1);


#clear anything that is currently in the analysis descriptions table

my $query= "delete from analysis_description";

$dbh->do($query);

$support->log("Deleted the current contents of the 'analysis_description' table.\n",1);

# fetch all the analysis-ids and logic names from the analysis table

$query= "select analysis_id, logic_name from analysis";


my $sth= $dbh->prepare($query);
$sth->execute;

my %name_mapper= ();


while(my($analysis_id, $logic_name)= $sth->fetchrow_array){

	$name_mapper{$logic_name}= $analysis_id;
	
}


$sth= $dbh->prepare(qq{insert into analysis_description (analysis_id, description, display_label) values (? , ?, '')});

foreach my $logic_name (keys %desc_mapper){

	my $analysis_id= $name_mapper{$logic_name};
	
	if(! defined($analysis_id)){
		
		print "Warning: No analysis-id for logic name: $logic_name\n";
		$support->log_warning("Unable to add description for logic_name: $logic_name, as no corresponding analysis_id in table analysis\n",2);
		
		next;
	
	
	}
	
	my $desc= $desc_mapper{$logic_name};
	$sth->execute($analysis_id, $desc);
	
	$support->log("Description added in table analysis_description for logic_name: $logic_name\n",1);
	
	
}




$support->finish_log;

