#!/usr/local/bin/perl

=head1 NAME

backup_vega_databases.pl - dump databases to a backed up location

=head1 SYNOPSIS

backup_vega_databases.pl [options]

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
    -i, --interactive=0|1               run script interactively (default: true)
    -n, --dry_run, --dry=0|1            don't write results to database
    -h, --help, -?                      print help (this message)

Specific options:

	  --dbfile=NAME                     name of file containing list database to backup
   	  --backup_directory=NAME           location of backed up files
	  --backup_host=NAME                location of backup directory

=head1 DESCRIPTION

Use this script to dump one or more databases using mysqldump and gzip to a location
that is backed up.

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
use Data::Dumper;

$| = 1;

our $support = new Bio::EnsEMBL::Utils::ConversionSupport($SERVERROOT);

# parse options
$support->parse_common_options(@_);
$support->parse_extra_options(
							  'dbfile=s',
                              'backup_directory=s',
                              'backup_host=s',
							 );
$support->allowed_params($support->get_common_params,
                         'dbfile',
                         'backup_directory',
                         'backup_host',
						);

if ($support->param('help') or $support->error) {
    warn $support->error if $support->error;
    pod2usage(1);
}

# ask user to confirm parameters to proceed
$support->confirm_params;

# get log filehandle and print heading and parameters to logfile
$support->init_log;

# connect to databases and get adaptors
my $dba = $support->get_database('ensembl');
my $dbh = $dba->dbc->db_handle;

#get databases to backup
my $tables = [];
if ($support->param('dbfile')) {
	my $IN = $support->filehandle('<', $support->param('dbfile'));
	while (<$IN>) {
		chomp;
		next if /^\s*$/;
		next if /^#/;
		push @$tables, [$_];
	}
}
else {
	my $sth = $dbh->prepare("show databases");
	$sth->execute;
	$tables = $sth->fetchall_arrayref;
}

my $host = $support->param('host');
my $port = $support->param('port');
my $user = $support->param('user');
my $pass = $support->param('pass');
my $backup_location = $support->param('backup_host').':'.$support->param('backup_directory');

#backup tables
foreach my $r (@$tables) {
	my $db = $r->[0];
	$support->log_stamped("\nBacking up $db...\n");
	my $filename = $db . '.sql.gz';
	my $r = `mysqldump --opt -u $user -p$pass -h $host -P $port $db | gzip > $filename`;
	$support->log("mysqldump: $r\n",1);
	$r = `scp $filename $backup_location`;
	$support->log("scp: $r\n",1);
	$r = `rm $filename`;
}

#warn Dumper($tables);
$support->finish_log;
