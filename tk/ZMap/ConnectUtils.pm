package ZMap::ConnectUtils;

use strict;
use warnings;
use Exporter;
use X11::XRemote;
#use XML::Simple;

our @ISA    = qw(Exporter);
our @EXPORT = qw(parse_params 
                 parse_response
                 $WAIT_VARIABLE);
our $WAIT_VARIABLE = 0;

sub parse_params{
    my ($pairs_string) = shift;
    my ($param, $value, $out, $tmp);

    $pairs_string =~ s/\s//g; # no spaces allowed.
    my (@pairs)   = split(/[;]/,$pairs_string);

    foreach (@pairs) {
	($param,$value) = split('=',$_,2);
	next unless defined $param;
	next unless defined $value;
	push (@{$tmp->{$param}},$value);
    }
    # this removes the arrays if there's only one element in them.
    # I'm not sure I like this behaviour, but it means it's not
    # necessary to remember to add ->[0] to everything. 
    $out = { 
        map { scalar(@{$tmp->{$_}}) > 1 ? ($_ => $tmp->{$_}) : ($_ => $tmp->{$_}->[0]) } keys(%$tmp) 
        };
    
    return $out;
}

sub parse_response{
    my $resposne = shift;
    my $delimit  = X11::XRemote::delimiter();

    my ($status, $xml) = split(/$delimit/, $resposne, 2);

    my $parser = '';#XML::Simple->new();
    my $hash   = $parser->XMLin($xml);
    
    return ($status, $hash);
}

1;
__END__

=pod

=head1 DESCRIPTION

 A bit like cgi parse_params

=cut
