package CanvasWindow::XRemoteExample;

=pod

=head1 NAME

CanvasWindow::XRemoteExample

=head1 DESCRIPTION

 Hopefully  a simple  example usage  of  ZMap::Connect.  ZMap::Connect
should make it easy to  integrate the connection method into any other
Tk application.  Famous last words I hear...

 Here I'm  using the  CanvasWindow module  as a base,  so a  number of
methods are found in that module which are used here.

=cut

use strict;
use warnings;
use base 'CanvasWindow';
use ZMap::Connect;
use Tk;


#=======================================================================
# Setup the object
#=======================================================================
sub new{
    my ($pkg, @args) = @_;
    my $self = $pkg->SUPER::new(@args);
    # Some setup for the User Interface....
    my $close_window = sub{
        $self->canvas->toplevel->destroy;
        $self = undef;
    };
    $self->canvas->toplevel->protocol('WM_DELETE_WINDOW',  $close_window);

    # Create a Button Frame for the buttons, etc...
    my $button_frame = $self->canvas->toplevel->Frame->pack(-side => 'top', -fill => 'x');
    # Add an Entry widget for the window ID
    my $label   = $button_frame->Label(-text => "Window ID:")->pack(-side => "left");
    my $idEntry = $button_frame->Entry(-width        => 20,
                                       -textvariable => $self->set_note_ref(),
                                       -font         => ['Helvetica', 12, 'normal'],
                                       )->pack(-side => 'left');
    # Setup some buttons to do things.
    my $zoomIn = $button_frame->Button(
                                       -text       => 'Zoom In',
                                       -command    => sub { 
                                           my $a = ${$self->set_note_ref};
                                           $self->current_window_id($a);
                                           $self->make_request('zoom_in');
                                       },
                                       )->pack(-side => 'right');
    my $zoomOut = $button_frame->Button(
                                        -text       => 'Zoom Out',
                                        -command    => sub { 
                                            my $a = ${$self->set_note_ref};
                                            $self->current_window_id($a);
                                            $self->make_request('zoom_out');
                                        },
                                        )->pack(-side => 'right');
    my $startTalking = $button_frame->Button(
                                        -text       => '"Connect"',
                                        -command    => sub { 
                                            my $a = ${$self->set_note_ref};
                                            $self->current_window_id($a);
                                            $self->getZMap2CreateClient;
                                        },
                                        )->pack(-side => 'right');


    # NOW THE INTERESTING BITS
    # Callback which gets called when the window is sent a message via
    # an atom.   The $c  and $r are  the ZMap::Connect object  and the
    # Request String respectively and  are passed by the ZMap::Connect
    # module to the callback.   Additional arguments maybe supplied as
    # data via a list ref as  the third argument to init (see below or
    # perldoc ZMap::Connect)
    my $callback = sub { 
        my ($c,$r,$canvas) = @_;
        my $rp = "txt was created on $canvas...";
        $canvas->createText(rand(100),rand(100),
                            -text   => "req: $r\nrep: $rp",
                            -anchor => 'nw');
        return (200,$rp);
    };
    # Create a new Object
    my $zmapConnector = ZMap::Connect->new(-server => 1);
    # Initialise it Args: Tk_Widget, callback, list ref of data for the callback.
    $zmapConnector->init($button_frame, $callback, [$self->canvas]);
    # store the new object for use later
    $self->zmap_connect($zmapConnector);

    # just a message.
    warn join(" ", $zmapConnector->server_window_id(), $zmapConnector->request_name, $zmapConnector->response_name) ;

    return $self;
}
#=======================================================================
# Implement something similar to these methods to keep hold of a couple 
# of objects/properties...
#=======================================================================
sub xclient{
    my ($self, $id) = @_;
    unless ($id){ warn "usage: $self->xclient(id);\n"; return; }
    my $client = $self->{'_xclients'}->{$id};
    if(!$client){
        # make a new client to send commands to a remote window
        $client = X11::XRemote->new(-id     => $id, 
                                    -server => 0,
                                    -_DEBUG => 1
                                    );
        # keep a hash of clients hashed on window id of remote window.
        $self->{'_xclients'}->{$id} = $client;
    }
    return $client;
}
sub zmap_connect{
    my ($self, $store) = @_;
    my $zmap = $self->{'_zmap'};
    if($store && !$zmap){
        my $id = $store->server_window_id();
        # keep the "client" with the others so we don't send ourself
        # messages. N.B. It is a SERVER, NOT a client
        $self->{'_xclients'}->{$id} = $store->xremote;
        # keep the object for later use.
        $self->{'_zmap'} = $store;
    }
    return $zmap;
}





#=======================================================================
# just some utility methods to make the module work, not essential to
# other implementations, but possibly useful.
#=======================================================================

sub getZMap2CreateClient{
    my ($self) = @_;
    my $id = $self->current_window_id();
    my $cnnct = $self->{'_connectedTo'}->{$id};
    if(!$cnnct){
        $self->xclient($id)->send_commands($self->zmap_connect->connect_request);
        $self->{'_connectedTo'}->{$id} = 1;
    }
}
sub set_note_ref{
    my ($self, $search) = @_;
    $self->{'_set_note'} = $search if $search;
    my $n = '';
    $self->{'_set_note'} ||= \$n;
    return $self->{'_set_note'};
}
sub current_window_id{
    my ($self, $id) = @_;
    $self->{'_cid'} = $id if $id;
    return $self->{'_cid'};
}
sub make_request{
    my ($self, @commands) = @_;
    my $xr = $self->xclient($self->current_window_id);
    push(@commands, 'zoom_in') unless @commands;
    my @a = $xr->send_commands(@commands);
    print "Sent commands\n";
    for(my $i = 0; $i < @commands; $i++){
        print $commands[$i] . " resulted in " . $a[$i] . "\n";
    }
}


1;
# remove next line to use bindDump()
__END__

# ====================================== #
# Just bindDump below.                   #
#  - I was using this 4 debugging comes  #
#    from Tk::Widget I think             #
# ====================================== #

sub bindDump {

    # Dump lots of good binding information.  This pretty-print subroutine
    # is, essentially, the following code in disguise:
    #
    # print "Binding information for $w\n";
    # foreach my $tag ($w->bindtags) {
    #     printf "\n Binding tag '$tag' has these bindings:\n";
    #     foreach my $binding ($w->bind($tag)) {
    #         printf "  $binding\n";
    #     }
    # }

    my ($w) = @_;

    my (@bindtags) = $w->bindtags;
    my $digits = length( scalar @bindtags );
    my ($spc1, $spc2) = ($digits + 33, $digits + 35);
    my $format1 = "%${digits}d.";
    my $format2 = ' ' x ($digits + 2);
    my $n = 0;

    my @out;
    push @out, sprintf( "\n## Binding information for '%s', %s ##", $w->PathName, $w );

    foreach my $tag (@bindtags) {
        my (@bindings) = $w->bind($tag);
        $n++;                   # count this bindtag

        if ($#bindings == -1) {
            push @out, sprintf( "\n$format1 Binding tag '$tag' has no bindings.\n", $n );
        } else {
            push @out, sprintf( "\n$format1 Binding tag '$tag' has these bindings:\n", $n );

            foreach my $binding ( @bindings ) {
                my $callback = $w->bind($tag, $binding);
                push @out, sprintf( "$format2%27s : %-40s\n", $binding, $callback );

                if ($callback =~ /SCALAR/) {
                    if (ref $$callback) {
                        push @out, sprintf( "%s %s\n", ' ' x $spc1, $$callback );
                    } else {
                        push @out, sprintf( "%s '%s'\n", ' ' x $spc1, $$callback );
                    }
                } elsif ($callback =~ /ARRAY/) {
                    if (ref $callback->[0]) {
                        push @out, sprintf( "%s %s\n", ' ' x $spc1, $callback->[0], "\n" );
                    } else {
                        push @out, sprintf( "%s '%s'\n", ' ' x $spc1, $callback->[0], "\n" );
                    }
                    foreach my $arg (@$callback[1 .. $#{@$callback}]) {
                        if (ref $arg) {
                            push @out, sprintf( "%s %-40s", ' ' x $spc2, $arg );
                        } else {
                            push @out, sprintf( "%s '%s'", ' ' x $spc2, $arg );
                        }
			
                        if (ref $arg eq 'Tk::Ev') {
                            if ($arg =~ /SCALAR/) {
                                push @out, sprintf( ": '$$arg'" );
                            } else {
                                push @out, sprintf( ": '%s'", join("' '", @$arg) );
                            }
                        }

                        push @out, sprintf( "\n" );
                    } # forend callback arguments
                } # ifend callback

            } # forend all bindings for one tag

        } # ifend have bindings

    } # forend all tags
    push @out, sprintf( "\n" );
    return @out;

} # end bindDump

1;

# ====================================== #
__END__
