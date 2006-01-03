#!/usr/local/bin/perl

=head1 NAME

gene_type_prefix_to_author.pl - convert gene.type prefixes to author suffixes

=head1 SYNOPSIS

gene_type_prefix_to_author.pl [options]

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

This script strips off prefixes from gene.type (or rather gene.biotype in the 
current schema) and assigns the affected genes to different authors. Author
names are created as the orininal author suffixed with a colon followed by the
gene.type prefix.

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
my $dba = $support->get_database('otter');
my $dbh = $dba->dbc->db_handle;

# get all gene.biotype prefix / author combinations
my $sql = qq(
    SELECT 
            g.biotype,
            a.author_id,
            a.author_name,
            a.author_email
    FROM
            gene g,
            gene_info gi,
            gene_stable_id gsi,
            author a
    WHERE   g.gene_id = gsi.gene_id
    AND     gsi.stable_id = gi.gene_stable_id
    AND     a.author_id = gi.author_id
    GROUP BY g.biotype, gi.author_id
);
my $sth1 = $dbh->prepare($sql);
$sth1->execute;

my %new_authors;
while (my ($type, $id, $name, $email) = $sth1->fetchrow_array) {
    $type =~ /(.*):(.*)/;

    # nothing to do if biotype doesn't contain a colon
    next unless $2;
    
    # create hash of distinct new author names
    $new_authors{"$name:$1"} = [$id, $email, $1];
}

$support->log("Looping over new authors...\n");
$sql = qq(INSERT INTO author values (?, ?, ?));
my $sth2 = $dbh->prepare($sql);
my $new_id = 2001;
foreach my $author (keys %new_authors) {
    my $id = $new_authors{$author}->[0];
    my $prefix = $new_authors{$author}->[2];

    # add new author
    $support->log("Adding new author (id, email, name) ".
        join(", ", $new_id, $new_authors{$author}->[1], $author)."\n", 1);
    unless ($support->param('dry_run')) {
        $sth2->execute($new_id, $new_authors{$author}->[1], $author);
    }
    $support->log("Done.\n", 1);
    
    # assign genes to new authors
    my $sql2 = qq(
        UPDATE
                gene g,
                gene_info gi,
                gene_stable_id gsi,
                author a
        SET 
                gi.author_id = $new_id
        WHERE   g.gene_id = gsi.gene_id
        AND     gsi.stable_id = gi.gene_stable_id
        AND     gi.author_id = a.author_id
        AND     a.author_id = $id
        AND     g.biotype like '$prefix:%'
    );
    $support->log("Assigning genes to new authors...\n", 1);
    if ($support->param('dry_run')) {
        $support->log_verbose("$sql2\n");
    } else {
        $dbh->do($sql2);
    }
    $support->log("Done.\n", 1);
    
    # assign transcripts to new authors
    my $sql3 = qq(
        UPDATE
                transcript t,
                transcript_info ti,
                transcript_stable_id tsi,
                author a
        SET 
                ti.author_id = $new_id
        WHERE   t.transcript_id = tsi.transcript_id
        AND     tsi.stable_id = ti.transcript_stable_id
        AND     ti.author_id = a.author_id
        AND     a.author_id = $id
        AND     t.biotype like '$prefix:%'
    );
    $support->log("Assigning transcripts to new authors...\n", 1);
    if ($support->param('dry_run')) {
        $support->log_verbose("$sql3\n");
    } else {
        $dbh->do($sql3);
    }
    $support->log("Done.\n\n", 1);

    $new_id++;
}
$support->log("Done.\n\n");

# strip gene.biotype and transcript.biotype prefix
$support->log("Stripping prefixes from gene.biotype and transcript.biotype...\n");
$sql = qq(
    UPDATE gene
    SET biotype = SUBSTRING_INDEX(biotype, ':', -1)
);
$dbh->do($sql) unless ($support->param('dry_run'));
$sql = qq(
    UPDATE transcript
    SET biotype = SUBSTRING_INDEX(biotype, ':', -1)
);
$dbh->do($sql) unless ($support->param('dry_run'));
$support->log("Done.\n\n");

# finish logfile
$support->finish_log;

