#!/usr/bin/perl

use strict;
use warnings;

use Digest::MD5 qw(md5_hex);
use URI::Escape qw(uri_escape uri_unescape);
use Time::HiRes qw(time);

my $gff_source;
sub _log_prefix {
    my $now = CORE::time;
    return sprintf "%02d:%02d:%02d: %6d: %-35s "
        , (localtime($now))[2,1,0], $$, ($gff_source || '???');
}

my $log_file;

sub log_message {
    my ($message) = @_;
    return unless $log_file;
    printf $log_file "%s: %s\n", _log_prefix, $message;
    return;
}

sub log_chunk {
    my ($prefix, $chunk) = @_;
    return unless $log_file;
    my $prefix_full = sprintf "%s: %s: ", _log_prefix, $prefix;
    chomp $chunk;
    $chunk .= "\n";
    $chunk =~ s/^/$prefix_full/gm;
    print $log_file $chunk;
    return;
}

my $LOG     = 1;

my %args;
foreach my $pair (@ARGV) {
    my ($key, $val) = split(/=/, $pair);
    $args{uri_unescape($key)} = uri_unescape($val);
}

$gff_source = $args{gff_source};

# test case
die "failing as required" if $args{'fail'};

# pull off arguments meant for us

my $url_root        = delete $args{'url_root'};
my $server_script   = delete $args{'server_script'};
my $session_dir     = delete $args{'session_dir'};
my $cookie_jar      = delete $args{'cookie_jar'};
my $process_gff     = delete $args{'process_gff_file'};

chdir($session_dir) or die "Could not chdir to '$session_dir'; $!";

open $log_file, '>>', 'gff_log.txt';
{
    ## no critic (InputOutput::ProhibitOneArgSelect)
    my $old_fh = select $log_file;
    $| = 1; ## no critic (Variables::RequireLocalizedPunctuationVars)
    select $old_fh;
}
log_message "starting";

$args{log} = 1 if $LOG; # enable logging on the server

# concatenate the rest of the arguments into a parameter string

my $params = join '&', map {
    uri_escape($_).'='.uri_escape($args{$_});
} sort keys %args;

my $gff_filename = sprintf '%s_%s.gff', $gff_source, md5_hex($params);

my $top_dir = 'gff_cache';

unless (-d $top_dir) {
    # Cannot check return value from mkdir() because another instance of the
    # script is likely to have made the directory since many run in parallel!
    mkdir $top_dir;
    unless (-d $top_dir) {
        die "Failed to create toplevel cache directory: $!\n";
    }
}

my $cache_file = $top_dir.'/'.$gff_filename;

if (-e $cache_file) {
    # cache hit
    log_message "cache file: $gff_filename: cache hit";
    open my $gff_file, '<', $cache_file or die "Failed to open cache file: $!\n";
    while (<$gff_file>) { print; }
    close $gff_file or die "Failed to close cache file: $!\n";
    close STDOUT or die "Error writing to STDOUT; $!";
    exit;
}

# cache miss

# only require these packages now, so we don't take the import hit on a cache hit

require LWP::UserAgent;
require HTTP::Request;
require HTTP::Cookies::Netscape;
require DBI;

log_message "cache file: $gff_filename: cache miss";

my $request = HTTP::Request->new;

$request->method('GET');

#$request->accept_decodable(HTTP::Message::Decodable);

my $url = $url_root . '/' . $server_script . '?' . $params;

log_message "http: URL: $url";

$request->uri($url);

# create a user agent to send the request

my $ua = LWP::UserAgent->new(
    timeout             => 9000,
    env_proxy           => 1,
    agent               => $0,
    cookie_jar          => HTTP::Cookies::Netscape->new(file => $cookie_jar),
    protocols_allowed   => [qw(http https)],
);

# do the request
log_message "http: start";
my $request_start_time = time;
my $response = $ua->request($request);
die "No response for $gff_source\n" unless $response;
my $request_finish_time = time;
my $request_time = time_diff($request_start_time, $request_finish_time);
log_message "http: finish: $request_time";

if ($response->is_success) {

    my $gff = $response->decoded_content;
    die "Unexpected response for $gff_source: $gff\n" unless $gff =~ /EnsEMBL2GFF/;

    # cache the result
    log_message "caching: start";
    my $cache_start_time = time;
    open my $cache_file_h, '>', $cache_file or die "Cannot write to cache file '$cache_file'; $!\n";
    print $cache_file_h $gff;
    close $cache_file_h or die "Error writing to '$cache_file'; $!";
    my $cache_finish_time = time;
    my $cache_time = time_diff($cache_start_time, $cache_finish_time);
    log_message "caching: finish: $cache_time";

    # update the SQLite db
    log_message "SQLite update: start";
    my $sqlite_update_start_time = time;
    my $dbh = DBI->connect("dbi:SQLite:dbname=$session_dir/otter.sqlite", undef, undef, {
        RaiseError => 1,
        AutoCommit => 1,
    });
    eval {
        $dbh->begin_work;
        my $sth = $dbh->prepare(
            q{ UPDATE otter_filter SET done = 1, failed = 0, gff_file = ?, process_gff = ? WHERE filter_name = ? }
            );
        $sth->execute($cache_file, $process_gff || 0, $gff_source);
    };
    if (my $err = $@) {
        $dbh->rollback;
        my $msg = "Update of otter_filter table in SQLite db failed; $err";
        log_message $msg;
        die $msg;
    }
    else {
        $dbh->commit;
    }
    $dbh->disconnect;
    my $sqlite_update_finish_time = time;
    my $sqlite_update_time = time_diff($sqlite_update_start_time, $sqlite_update_finish_time);
    log_message "SQLite update: finish: $sqlite_update_time";

    # if (rand() < 0.5) {
    #     die "Horribly";
    # }
    
    # Send data to zmap on STDOUT
    log_message "sending data: start";
    my $send_data_start_time = time;
    print STDOUT $gff;
    my $send_data_finish_time = time;
    my $send_data_time = time_diff($send_data_start_time, $send_data_finish_time);
    log_message "sending data: finish: $send_data_time";

    # zmap waits for STDOUT to be closed as an indication that all
    # data has been sent, so we close the handle now so that zmap
    # doesn't tell otterlace about the successful loading of the column
    # before we have the SQLite db updated and the cache file saved.
    close STDOUT or die "Error writing to STDOUT; $!";
}
else {

    my $res = $response->content;
    log_chunk 'http: error', $res;

    my $err_msg;

    if (my ($err) = $res =~ /ERROR:[[:space:]]*(.+)/s) {
        $err =~ s/\A(^-+[[:blank:]]*EXCEPTION[[:blank:]]*-+\n)+//m; # remove boring initial lines
        $err =~ s/\n.*//s; # keep only the first line
        $err_msg = $err;
    }
    elsif ($res =~ /The Sanger Institute Web service you requested is temporarily unavailable/) {
        my $code = $response->code;
        my $message = $response->message;
        $err_msg = "This Sanger web service is temporarily unavailable: status = ${code} ${message}";
    }
    else {
        $err_msg = $res;
    }

    die "Webserver error for $gff_source: $err_msg\n";
}

sub time_diff {
    my ($start_time, $end_time) = @_;
    
    return sprintf "time (sec): %.3f", $end_time - $start_time;
}
