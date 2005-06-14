package ZMap::ConnectUtils;

use strict;
use warnings;
use Exporter;
use X11::XRemote;
use XML::Simple;


our @ISA    = qw(Exporter);
our @EXPORT = qw(parse_params 
                 parse_response
                 $WAIT_VARIABLE);
our @EXPORT_OK = qw(xclient_with_id
                    xclient_with_name
                    list_xclient_names
                    delete_xclient_with_id
                    delete_xclient_with_name);
our %EXPORT_TAGS = ('caching' => [@EXPORT_OK],
                    'all'     => [@EXPORT, @EXPORT_OK]
                    );
our $WAIT_VARIABLE = 0;

=pod

=head1 NAME

ZMap::ConnectUtils

=head1 DESCRIPTION

 functions which might get written over and over again.

=head1 PARSING FUNCTIONS

=head2 parse_params(string)

parse a string like name = value ; name = valU ; called = valid
into 
 {
    name   => [qw(value valU)],
    called  => 'valid'
 }

=cut

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

=head2 parse_response(string)

parse a response into (status, hash of xml)

=cut

sub parse_response{
    my $response = shift;
    my $delimit  = X11::XRemote::delimiter();

    my ($status, $xml) = split(/$delimit/, $response, 2);

    my $parser = XML::Simple->new();
    my $hash   = $parser->XMLin($xml);
    
    return wantarray ? ($status, $hash) : $hash;
}


{
    # functions to cache clients
    my $CACHED_CLIENTS = {
        #'0x000000' => ['<X11::XRemote>', 'name']
    };
    my $AUTO_INCREMENT = 0;

=head1 FUNCTIONS TO CACHE CLIENTS

=over 5

=item I<create and retrieve>

=back

=head2 xclient_with_id(id)

return the xclient with specified id, creating if necessary.

=cut

    sub xclient_with_id{
        my ($id) = @_;
        # user can't use a name for some reason, lets make one up
        my $name = __PACKAGE__ . $AUTO_INCREMENT++;
        return xclient_with_name($name, $id);
    }

=head2 xclient_with_name(name, [id])

return the xclient with specified name, creating if id supplied.

=cut

    sub xclient_with_name{
        my ($name, $id) = @_;
        if(!$id){
            # we need to look it up
            ($id) = grep { $CACHED_CLIENTS->{$_}->[1] eq $name } keys %$CACHED_CLIENTS;
        }else{ 
            # check we haven't already got that name
            if(my ($cachedid) = grep { $CACHED_CLIENTS->{$_}->[1] eq $name} keys %$CACHED_CLIENTS){
                unless($id eq $cachedid){
                    warn "we've already got name $name with id $cachedid".
                        " cannot add id $id, try delete_xclient_with_id($id) first\n";
                    return undef;
                }
            }
        }
        return unless $id;

        local *xclient = sub{
            my ($id, $name)   = @_;
            my $client = $CACHED_CLIENTS->{$id}->[0];
            if(!$client){
                $client = X11::XRemote->new(-id     => $id, 
                                            -server => 0,
                                            -_DEBUG => 1
                                            );
                $CACHED_CLIENTS->{$id} = [ $client, $name ];
            }

            return $client;
        };

        return xclient($id, $name);
    }

    sub list_xclient_names{
        flush_bad_windows();
        my @list = ();
        foreach my $obj_name(values(%$CACHED_CLIENTS)){
            next if $obj_name->[0]->_is_server();
            push(@list, $obj_name->[1]);
        }
        return @list;
    }

=over 5

=item I<removal>

=head2 delete_xclient_with_id(id)

remove the xclient with specified id.

=cut

    # remove
    sub delete_xclient_with_id{
        my ($id) = @_;
        delete $CACHED_CLIENTS->{$id};
    }

=head2 delete_xclient_with_name(name)

remove the xclient with specified name.

=cut

    sub delete_xclient_with_name{
        my ($name) = @_;
        my ($id) = grep { $CACHED_CLIENTS->{$_}->[1] eq $name } keys %$CACHED_CLIENTS;
        delete_xclient_with_id($id);
    }

    sub flush_bad_windows{
        foreach my $id(keys %$CACHED_CLIENTS){
            my $obj = xclient_with_id($id);
            $obj->ping || delete_xclient_with_id($id);
        }
    }
}


1;
__END__

