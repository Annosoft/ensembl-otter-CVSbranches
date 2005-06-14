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
use ZMap::Connect qw(:all);
use Tie::Watch;
use Tk;
use Tk::BrowseEntry;

my $WATCH_METHOD = 0;
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
    my $kill = $button_frame->Button(-text => 'End Test',
                                     -command => $close_window
                                     )->pack(-side => 'left');
    my $idEntry = $button_frame->BrowseEntry(-label    => "Window ID:",
                                             -width    => 20,
                                             -variable => $self->entry_ref(),
                                             -font     => ['Helvetica', 12, 'normal'],
                                             -listcmd  => sub {
                                                 my $w = shift;
                                                 $w->choices([ list_xclient_names() ])
                                                 },
                                             )->pack(-side => 'left');

    $self->canvas->Tk::bind('<Control-Button-1>', 
                        sub {$self->make_request("newZmap seq = b0250 ; start = 1 ; end = 0")}) if !$WATCH_METHOD;

    # Setup some buttons to do things.
    my $open = 
        $button_frame->Button(-text       => 'open',
                              -command    => sub { 
                                  $self->make_request("newZmap seq = b0250 ; start = 1 ; end = 0");
                              })->pack(-side => 'right');
    my $zoomIn = 
        $button_frame->Button(-text       => 'Zoom In',
                              -command    => sub { 
                                  $self->make_request('zoom_in');
                              })->pack(-side => 'right');
    my $zoomOut = 
        $button_frame->Button(-text       => 'Zoom Out',
                              -command    => sub { 
                                  $self->make_request('zoom_out');
                              })->pack(-side => 'right');
    my $startTalking = 
        $button_frame->Button(-text       => 'start ZMap',
                              -command    => sub { 
                                  $self->start_zmap();
                              })->pack(-side => 'right');
    # Create a new Object
    my $zmapConnector = ZMap::Connect->new(-server => 1);
    # Initialise it Args: Tk_Widget, callback, list ref of data for the callback.
    $zmapConnector->init($button_frame, \&RECEIVE_FILTER, [$self, qw(simpleton register_client simple default ping)]);
    # store the new object for use later
    $self->zmapConnectObj($zmapConnector);

    # just a message.
    print join(" ", $zmapConnector->server_window_id(), $zmapConnector->request_name, $zmapConnector->response_name) . "\n" ;

    return $self;
}

#=====================================================================
# Button Functions
#=====================================================================
sub start_zmap{
    my ($self) = @_;
    return if xclient_with_name('main');
    if(my $pid = fork()){
        print "forked ok\n";
        return;
    }else{
        my @exe = qw(zmap --win_id);
        close STDERR if $^O eq 'linux';
        exec(@exe, $self->zmapConnectObj()->server_window_id());
    }
}
sub make_request{
    my ($self, @commands) = @_;

    my $xr = $self->current_xclient;
    unless($xr){
        warn "No current window.";
        return;
    }

    push(@commands, 'zoom_in') unless @commands;
    my @a = $xr->send_commands(@commands);
    my $canvasMessage = "Sent commands\n";
    for(my $i = 0; $i < @commands; $i++){
        $canvasMessage .= " - cmd: '" . $commands[$i] . "',  result: '" . $a[$i] . "'\n";
        my ($status, $n) = parse_response($a[$i]);
        my ($name, $id) = ($n->{zmapid}, $n->{windowid});
        if($status == 412){
            delete_xclient_with_id($n->{'error'}->{'windowid'});
        }
        if($name){
            xclient_with_name($name, $id) if $id;
            $self->set_entry_value($name);
        }
    }
    $self->write_on_canvas($canvasMessage);
}

#=======================================================================
# holds the ZMap::Connect object
#=======================================================================
sub zmapConnectObj{
    my ($self, $store) = @_;
    my $zmap = $self->{'_zmap'};
    if($store && !$zmap){
        # keep the object for later use.
        $self->{'_zmap'} = $store;
    }
    return $zmap;
}
# and sends a register_client
sub getZMap2CreateClient{
    my ($self) = @_;
    my $id = $self->current_xclient()->window_id();
    my $cnnct = $self->{'_connectedTo'}->{$id};
    if(!$cnnct){
        my $req = $self->zmapConnectObj->connect_request;
        $self->make_request($req);
        $self->{'_connectedTo'}->{$id} = 1;
    }
}

#=====================================================================
# RECEIVE_FILTER && functions it calls
#=====================================================================
sub RECEIVE_FILTER{ 
    my ($_connect, $_request, $_obj, @valid_list) = @_; 
    my ($_status, $_response) = (404, "<error>unknown command</error>");

    foreach my $valid(@valid_list){
        if($_request =~ s/^$valid//){
            ($_status, $_response)
                = $_obj->$valid($_request) if $_obj->can($valid);
            last;
        }
    }

    $_obj->write_on_canvas("req: $_request \nresp: $_response");

    return ($_status, $_response) ;    

}

sub register_client{
    my ($self, $request) = @_;

    my $p = parse_params($request);

    unless($p->{'id'} 
           && $p->{'request'}
           && $p->{'response'}){
        warn "mismatched request for register_client:\n", 
        "id, request and response required\n",
        "Got '$request'\n";
        return (333,333);
    }
    # now we're all ok we can do the registering
    my $client = xclient_with_name('main', $p->{'id'});
    $client->request_name( $p->{'request'} );
    $client->response_name($p->{'response'});

    if($WATCH_METHOD){
        my $watch = Tie::Watch->new(-variable => \$WAIT_VARIABLE,
                                    -debug    => 0,
                                    -store    => [ \&watching_sub, $self ],
                                    );
    }else{
        $self->canvas->eventGenerate('<Control-Button-1>', -when => 'tail');
    }
    $self->set_entry_value('main');

    return (200, "<message>accepted</message>");
}
sub watching_sub{
    my ($watch, $value) = @_;
    my ($self) = @{$watch->Args('-store')};
    my $xr     = xclient_with_name('main');
    my @clones = qw(b0250);

    $self->make_request("newZmap seq = @clones ; start = 1 ; end = 0");

    $watch->Unwatch;
}



#=====================================================================
# short helper bits
#=====================================================================
sub current_xclient{
    my ($self) = @_;
    return xclient_with_name(${$self->entry_ref});
}
sub entry_ref{
    my ($self) = @_;
    my $n = '';
    $self->{'_entry_ref'} ||= \$n;
    return $self->{'_entry_ref'};
}
sub set_entry_value{
    my ($self, $value) = @_;
    my $ref = $self->entry_ref();
    $$ref   = $value;
}
sub write_on_canvas{
    my ($self, $message) = @_;
    $self->canvas->delete('all');
    $self->canvas->createText(10, 10,
                              -text   => "$message",
                              -anchor => 'nw');

}

#=====================================================================
# Just to check for memory leaks
#=====================================================================
sub DESTROY{
    my $self = shift;
    warn "Destroying $self";
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
