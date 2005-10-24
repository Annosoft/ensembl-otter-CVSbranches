#!/usr/local/bin/perl

=head1 NAME

add_assembly.pl

=head1 SYNOPSIS

add_assembly.pl [options]

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
    --ensembldbname=NAME                use Ensembl database NAME
    --ensemblhost=HOST                  use Ensembl database host HOST
    --ensemblport=PORT                  use Ensembl database port PORT
    --ensembluser=USER                  use Ensembl database username USER
    --ensemblpass=PASS                  use Ensembl database passwort PASS

=head1 DESCRIPTION


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
    unshift(@INC, "$SERVERROOT/ensembl-otter/modules");
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
$support->parse_extra_options(
    'ensemblhost=s',
    'ensemblport=s',
    'ensembluser=s',
    'ensemblpass=s',
    'ensembldbname=s',
);
$support->allowed_params(
    $support->get_common_params,
    'ensemblhost',
    'ensemblport',
    'ensembluser',
    'ensemblpass',
    'ensembldbname',
);

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

my $s_host    = 'ecs2b';
my $s_user    = 'ensro';
my $s_pass    = '';
my $s_dbname  = 'ens_NCBI_31';
my $s_port    = 3306;

my $t_host    = 'ecs1d';
my $t_user    = 'ensadmin';
my $t_pass    = 'ensembl';
my $t_port    = 19322;
my $t_dbname  = 'otter_chr22_with_anal';

my $chr      = '22';
my $chrstart = 1;
my $chrend   = 100000000;
my $path     = 'NCBI31';

my $no_offset;

&GetOptions( 's_host:s'    => \$s_host,
             's_user:s'    => \$s_user,
             's_pass:s'    => \$s_pass,
             's_port:s'    => \$s_port,
             's_dbname:s'  => \$s_dbname,
             't_host:s'=> \$t_host,
             't_user:s'=> \$t_user,
             't_pass:s'=> \$t_pass,
             't_port:s'=> \$t_port,
             't_dbname:s'  => \$t_dbname,
             'chr:s'     => \$chr,
             'chrstart:n'=> \$chrstart,
             'chrend:n'  => \$chrend,
             'path:s'  => \$path,
	     'no_offset' => \$no_offset,
            );



my $sdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $s_host,
					       -user => $s_user,
					       -pass => $s_pass,
					       -port => $s_port,
					       -dbname => $s_dbname);

$sdb->assembly_type($path);

my $tdb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $t_host,
					       -user => $t_user,
					       -pass => $t_pass,
					       -port => $t_port,
					       -dbname => $t_dbname);

my $chr_sth = $tdb->prepare( qq {
   select chromosome_id from chromosome 
   where chromosome.name='$chr'
});

$chr_sth->execute;

my $chr_hashref = $chr_sth->fetchrow_hashref;
if (!defined $chr_hashref) {
  die "Couldn't find chromosome $chr in target database\n";
}

my $t_chrid = $chr_hashref->{'chromosome_id'};

print "Target chromosome id = $t_chrid\n";

# First check for duplicate clones used in contigs
my $dup_sth = $sdb->prepare( qq {
  select contig.name from 
         assembly, contig, dna, chromosome 
  where  assembly.chromosome_id=chromosome.chromosome_id and
         chromosome.name='$chr' and 
         contig.contig_id=assembly.contig_id and 
         assembly.chr_start>= $chrstart and
         assembly.chr_end <=$chrend and
         assembly.type = "$path" and
         contig.dna_id=dna.dna_id
  order by chr_start
});
        
$dup_sth->execute;


my %clones;
my %dupclones;
my $ndup = 0;
my $hashref;
while ($hashref = $dup_sth->fetchrow_hashref()) {
  my $contigname = $hashref->{'name'};
  my @bits = split /\./,$contigname;

  my $clone = $bits[0] . "." . $bits[1];

  if (exists($clones{$clone})) {
    print "Duplicate clone "  . $clone . "\n";
    $dupclones{$clone} = 1;
    $ndup++;
  }

  $clones{$clone} = $_;
}



my $sth = $sdb->prepare(qq {
  select chromosome.chromosome_id,
         chr_start,
         chr_end,
         superctg_name,
         superctg_start,
         superctg_end,
         superctg_ori,
         contig.name,
         contig.clone_id,
         contig_start,
         contig_end,
         contig_ori,
         type 
  from assembly,contig,chromosome 
  where assembly.type = "$path" and
        assembly.chromosome_id=chromosome.chromosome_id and
        chromosome.name = '$chr' and
        assembly.chr_start>= $chrstart and
        assembly.chr_end <=$chrend and
        assembly.contig_id = contig.contig_id
  order by chr_start
});
        
$sth->execute;
my $rca = $tdb->get_RawContigAdaptor();

my $t_ca = $tdb->get_CloneAdaptor();
my $s_ca = $sdb->get_CloneAdaptor();

my $hashref;
my $ndiffseq = 0;
my $totfailed = 0;
my $nnotfound = 0;
my $nwritten = 0;
my $ndupcontig = 0;
while ($hashref = $sth->fetchrow_hashref()) {

  my $contigname = $hashref->{'name'};
  my @bits = split /\./,$contigname;

  my $clonename = $bits[0] . "." . $bits[1];

  $hashref->{'t_contig_start'}=$hashref->{'contig_start'};
  $hashref->{'t_contig_end'}=$hashref->{'contig_end'};

  if (!exists($dupclones{$clonename})) {

    my $contig = $rca->fetch_by_name($hashref->{'name'});
   
    # check for contigs with name without version (AC.V.ST.ED -> AC.ST.ED)
    if (!defined($contig)) {
      print STDOUT "Missing contig " . $hashref->{'name'} . " - trying without version\n";
      my $acc = $hashref->{'name'};
      $acc =~ s/(\S+)\.[0-9]+(\.[0-9]+\.[0-9]+)$/\1\2/;
      #print $hashref->{'name'}." ".$acc."\n";
      #exit 0;
      $contig = $rca->fetch_by_name($acc);
      if(defined($contig)){
	print "Looking for $acc successful\n";
      }
    }

    # check for contigs with without start-end (AC.V.ST.ED -> AC.V)
    if (!defined($contig)) {
      print STDOUT "Missing contig " . $hashref->{'name'} . " - trying clone\n";
      my $acc = $hashref->{'name'};
      $acc =~ s/\.[0-9]*\.[0-9]*$//;
      my $sth = $tdb->prepare("select name from contig where contig.name like '$acc%'");
      $sth->execute;
      my $hashref2 = $sth->fetchrow_hashref;
      if (defined($hashref2)) {
	my $name2=$hashref2->{'name'};
        print "Looking for $acc did find $name2\n";
	$contig = $rca->fetch_by_name($hashref2->{'name'});

	# contig fetch worked - now fix start end
	my $name = $hashref->{'name'};
	if($name=~/^(\w+\.\d+)\.(\d+)\.(\d+)$/ && !$no_offset){
	  # change from AC.V.ST.ED -> AC.V.1.ED
	  my($sv,$st,$ed)=($1,$2,$3);
	  my $cst = $hashref->{'contig_start'};
	  my $ced = $hashref->{'contig_end'};
	  if($cst==1 && $ced==$ed-$st+1){
	    $cst=$st;
	    $ced=$ced+$st-1;
	    print STDOUT "  assembly coordinates remapped to $cst-$ced\n";
	    $hashref->{'t_contig_start'}=$cst;
	    $hashref->{'t_contig_end'}=$ced;
	  }
	}
      }
    }

    if (!defined($contig)) {
  
      # Fetch clone from source db
#      my $clone = $s_ca->fetch_by_dbID($hashref->{'clone_id'});
#  
#      foreach my $s_contig (@{$clone->get_all_Contigs}) {
#        my $seq = $s_contig->seq;
#        $s_contig->adaptor($rca);
#        $s_contig->seq($seq);
#      }
#  
#      # Write into target db
#      $clone->adaptor($t_ca);
#  
#      # Now store the clone
#      $t_ca->store($clone);
#  
#      $contig = $rca->fetch_by_name($hashref->{'name'});
    }
  
    if (defined($contig)) {
      my $s_clone = $s_ca->fetch_by_dbID($hashref->{'clone_id'});
      my $s_contig = @{$s_clone->get_all_Contigs}[0];

      my $s_assseq = $s_contig->subseq($hashref->{contig_start},$hashref->{contig_end});
      my $t_assseq = $contig->subseq($hashref->{t_contig_start},$hashref->{t_contig_end});
      if ($s_assseq ne $t_assseq) {
	print "start = " . $hashref->{contig_start} . " end " . $hashref->{contig_end} . "\n";
        print "NOTE Contig subseqs sequences different for $clonename\n";
        compare_seqs($s_assseq, $t_assseq);
        $ndiffseq++; 
        $totfailed++;
      } else {

        my $contig_id = $contig->dbID;
        print "Found contig " . $contig->name . "\n";
        my $sth2 = $tdb->prepare( qq:
          insert into assembly(chromosome_id,
                               chr_start,
                               chr_end,
                               superctg_name,
                               superctg_start,
                               superctg_end,
                               superctg_ori,
                               contig_id,
                               contig_start,
                               contig_end,
                               contig_ori,
                               type) 
                       values( 
                           $t_chrid,
                           $hashref->{chr_start},
                           $hashref->{chr_end},
                           "$hashref->{superctg_name}",
                           $hashref->{superctg_start},
                           $hashref->{superctg_end},
                           $hashref->{superctg_ori},
                           $contig_id,
                           $hashref->{t_contig_start},
                           $hashref->{t_contig_end},
                           $hashref->{contig_ori},
                           "$hashref->{type}"
                         ):);
        #print $sth2->{Statement} . "\n";
        $sth2->execute;
        $nwritten++;
      }
    } else {
      print "Didn't find " . $hashref->{'name'} . "\n";
      $nnotfound++; 
      $totfailed++;
    } 
  } else {
    print "Clone with multiple contigs $clonename\n";
    $totfailed++;
    $ndupcontig++;
  }
}

print "Number of assembly elements written = $nwritten\n";
print "Total failures = $totfailed\n";
print "No. with different seq = $ndiffseq\n";
print "No. not found in target db = $nnotfound\n";
print "No. of clones with duplicates = " . scalar(keys(%dupclones)) . " (n contig  = $ndupcontig)\n";

# finish logfile
$support->finish_log;


sub compare_seqs {
  my ($seq1, $seq2) = @_;

  $seq1 =~ s/(.{80})/$1\n/g;
  $seq2 =~ s/(.{80})/$1\n/g;

  # print "Chr = $chrstr\n";
  # print "Contig = " . $contigsubstr . "\n";

  if ($seq1 ne $seq2) {
    my $ndiffline = 0;
    my @seq2lines = split /\n/,$seq2;
    my @seq1lines = split /\n/,$seq1;
    for (my $linenum = 0; $linenum<scalar(@seq2lines); $linenum++) {
      if ($seq1lines[$linenum] ne $seq2lines[$linenum]) {
        $ndiffline++;
      }
    }
    print "N diff line = $ndiffline N line = " . scalar(@seq2lines)."\n";
    if ($ndiffline > 0.95*scalar(@seq2lines)) {
#        print "Chr = $chrstr\n";
#        print "Contig = " . $seq1 . "\n";
      for (my $linenum = 0; $linenum<scalar(@seq2lines); $linenum++) {
        if ($seq1lines[$linenum] eq $seq2lines[$linenum]) {
          print "Matched line in very different: $seq1lines[$linenum]\n";
        }
      }
    }
  }
}
