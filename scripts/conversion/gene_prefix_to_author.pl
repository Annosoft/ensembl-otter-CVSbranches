#!/usr/local/bin/perl

=head1 NAME

gene_type_prefix_to_author.pl - generate author_ids from gene.type prefixes

=head1 SYNOPSIS

human_gene_prefix_to_author.pl [options]

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

=head1 DESCRIPTION

This script uses the prefix from the gene display_xref to assign the correct author_id
to gene_info and transcript_info. It also updates author and author_group according to 
the %authors hash. As such it must be run after add_vega_xrefs and before otter_to_ensembl.pl

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Steve Trevanion <st3@sanger.ac.uk>

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
$support->allowed_params($support->get_common_params);

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
my $sa = $dba->get_SliceAdaptor();
my $ga = $dba->get_GeneAdaptor();

#details of authors (gene_prefix => [name, email, author_id]_ 
my %authors = (
			BCM    => ['BCM','sabo@bcm.tmc.edu',1 ],
			C22    => ['Sanger','chr22@sanger.ac.uk',2],
			GC     => ['Genoscope','ralph@genoscope.cns.fr',3 ],
			GD     => ['EGAG','egag@sanger.ac.uk',4 ],
			JGI    => ['JGI','UHellsten@lbl.gov',5 ],
            MIT    => ['Broad','daved@broad.mit.edu',6 ],
			SK     => ['Sick_Kids','jmacdonald@sickkids.ca',7 ],
			WU     => ['Washu','jspieth@watson.wust',8 ],
			HAVANA => ['Havana','vega@sanger.ac.uk',9],
			KO     => ['Havana','vega@sanger.ac.uk',9],
);


#prepare sql
my $sql = $dbh->prepare(qq(
        UPDATE  gene_info
        SET     author_id = ?
        WHERE   gene_stable_id = ?
));

my $sql2 = $dbh->prepare(qq(
          INSERT into author 
          VALUES (?,?,?,?)
));
my $sql3 = $dbh->prepare(qq(
          INSERT into author_group
          VALUES (?,?,?)
));
my $sql4 = $dbh->prepare(qq(
           UPDATE  transcript_info
           SET     author_id = ?
           WHERE   transcript_stable_id = ?
));

if (! $support->param('dry_run')) {
	##replace existing author and author_group tables with new stuff
	$support->log("Replacing data in author and author_group tables\n"); 		
	$dbh->do("delete from author_group");
	$dbh->do("delete from author");
	foreach my $k (keys %authors) {
		next if ($k eq 'KO');
		my $a = $authors{$k};
		$sql2->execute($a->[2],$a->[1],$a->[0],$a->[2]);
		$sql3->execute($a->[2],$a->[1],$a->[0]);
	}
}

foreach my $chr ($support->sort_chromosomes) {
    $support->log_stamped("> Chr $chr\n");
    my $slice = $sa->fetch_by_region('chromosome', $chr);
    foreach my $gene (@{ $ga->fetch_all_by_Slice($slice) }) {
		my $gid = $gene->stable_id;
		my $g_dbID = $gene->dbID;
		my $name = $gene->display_xref->display_id;
		my $det;			
		my $prefix = $name;
		#set authors for prefixed genes
		if ($prefix =~ s/(.*?):.*/$1/){
			if ($det = $authors{$prefix}) {
				$support->log("Updating gene $gid ($name) with author ".$det->[0]."\n"); 		
				if (! $support->param('dry_run')) {				
					$sql->execute($det->[2],$gid);
				}
			}
			else {
				$det = $authors{'HAVANA'};
				$support->log_warning("Setting gene $gid ($name) to default author ".$det->[0]."\n");
     			if (! $support->param('dry_run')) {					
			    	$sql->execute($det->[2],$gid);
				}
			}
		}
		else {
			#set Havana as author for non prefixed genes
			$det = $authors{'HAVANA'};
			$support->log("Setting gene $gid ($name) to default author ".$det->[0]."\n");
			if (! $support->param('dry_run')) {		
				$sql->execute($det->[2],$gid);
			}
		}

		foreach my $trans (@{$gene->get_all_Transcripts}) {
			my $tid = $trans->stable_id;
			$support->log_verbose("Setting author for transcript $tid to the same as the parent gene: ".$det->[0]."\n",1);
			if (! $support->param('dry_run')) {				
				$sql4->execute($det->[2],$tid);
			}
		}	
	}
}

$support->log("Done.\n\n");

# finish logfile
$support->finish_log;

