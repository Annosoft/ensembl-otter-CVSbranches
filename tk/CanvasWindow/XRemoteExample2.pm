package CanvasWindow::XRemoteExample2;

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
use Carp;
use Data::Dumper;
use Sys::Hostname;
use Hum::Ace::LocalServer;

my $WATCH_METHOD = 0;
#=======================================================================
# Setup the object
#=======================================================================
sub new{
    my ($pkg, @args) = @_;
    my $self = $pkg->SUPER::new(@args);
    # Some setup for the User Interface....
    my $close_window = sub{
        my $top = $self->canvas->toplevel;
        $top->destroy;
        reset_sigCHLD();
        $self = undef;
    };
    $self->canvas->toplevel->protocol('WM_DELETE_WINDOW',  $close_window);
    # Create a Button Frame for the buttons, etc...
    my $button_frame = $self->canvas->toplevel->Frame->pack(-side => 'top', -fill => 'x');
    # Add an Entry widget for the window ID


    my $kill = $button_frame->Button(-text => 'End Test',
                                     -command => sub {
                                         my $top = $self->canvas->toplevel;
                                         $top->withdraw;
                                     }
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
                        sub {$self->make_request("newZmap seq = 1.141003345-141198932 ; start = 1 ; end = 0")}) if !$WATCH_METHOD;
    $self->canvas->Tk::bind('<Destroy>', sub {$self = undef});
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
    my $query = sprintf('feature_find type = %s ; method = %s ; start = %d ; end = %d ; strand = %s ; feature = %s',
                        qw(transcript Processed_pseudogene 124583 125315 - RP11-417J8.4-001));
#                           RP11-417J8.4-001   -   124583   125315   TRANSCRIPT   Processed_pseudogene
    my $find = 
        $button_frame->Button(-text       => 'find RP11-417J8.4-001',
                              -command    => sub { 
                                  $self->make_request($query);
                              })->pack(-side => 'right');
    my $register = 
        $button_frame->Button(-text       => 'Register',
                              -command    => sub { 
                                  $self->getZMap2CreateClient();
                              })->pack(-side => 'right');

    my $startTalking2 = 
        $button_frame->Button(-text       => 'start sgiface & zmap',
                              )->pack(-side => 'right');
    $startTalking2->bind('<ButtonRelease-1>', [ sub { $self->start_sgiface_zmap } ]);
    $startTalking2->bind('<Destroy>', sub {$self = undef});

    # Create a new Object
    my $zmapConnector = ZMap::Connect->new(-server => 1);
    # Initialise it Args: Tk_Widget, callback, list ref of data for the callback.
    $zmapConnector->init($button_frame, \&RECEIVE_FILTER, [$self, qw(simpleton register_client simple default ping)]);
    # store the new object for use later
    $self->zmapConnectObj($zmapConnector);

    # just a message.
    print join(" ", $zmapConnector->server_window_id(), $zmapConnector->request_name, $zmapConnector->response_name, sprintf("0x%lx",$self->canvas->toplevel->wrapper)) . "\n" ;

    return $self;
}

#=====================================================================
# Button Functions
#=====================================================================

sub start_zmap{
    my ($self, @extra) = @_;

    if(xclient_with_name('main')){
        warn "zmap already started";
        return 0;
    }
    my $zmap    = 'zmap';
    my $wid     = $self->zmapConnectObj()->server_window_id();
    my @command = ($zmap, '--win_id', $wid, @extra);

    warn "Looking to fork and exec cmd @command \n";

    my $cb_wipe = sub {
        warn "pid: @_\n";
        flush_bad_windows();
        $self->write_on_canvas( Dumper $self->{'_children'} );
    };

    fork_exec(\@command, \%{$self->{'_children'}}, 0, $cb_wipe) || warn "failed";
}

sub start_sgiface_zmap{
    my ($self) = @_;
    if(xclient_with_name('main')){
        warn "already started";
        return 0;
    }
    if($self->spawn_sgifaceserver()){
        print "Good I can connect to sgifaceserver\n";
        $self->write_dot_zmap();
        $self->start_zmap('--conf_dir', $self->{'_SGIF_LOCAL_SERVER'}->path . '/.ZMap');
    }else{
        warn "BooHoo I can't connect\n";
    }
}
sub spawn_sgifaceserver{
    my ($self) = @_;
    my $sgif   = $self->{'_SGIF_LOCAL_SERVER'};
    my $path   = '/nfs/team71/analysis/rds/demo/test-data/lace.775074.ro_1';
    unless($sgif){
        $sgif = Hum::Ace::LocalServer->new($path);
#        $sgif->server_executable("/nfs/team71/acedb/edgrif/acedb/CODE/acedb/bin.ALPHA_5/sgifaceserver");
        $sgif->server_executable("sgifaceserver");
        $sgif->timeout_string('0:0');
        $sgif->start_server() or return 0; # this only check the fork was successful
        $sgif->ace_handle(1)  or return 0; # this checks it can connect
        $self->{'_SGIF_LOCAL_SERVER'} = $sgif;
    }
    return 1;
}
sub write_dot_zmap{
    my ($self) = @_;
    my $dir    = $self->{'_SGIF_LOCAL_SERVER'}->path;
    my $file   = "${dir}/.ZMap/ZMap";
    my $fh;
    eval{
        # directory should be made already
        open($fh, ">$file") or die "write_dot_zmap: error writing file '$file', $!";
    };
    warn "Error in :$@" if $@;
    unless($@){
        my $content = $self->dot_zmap_content();
        print $fh $content;
        return 1;
    }
    return 0;
}
sub dot_zmap_content{
    my ($self) = @_;
    my $fmt = qq`source\n{\nurl = "%s://%s:%s@%s:%d"\nsequence = %s\nwriteback = %s\n#stylesfile = "%s"\n}\n`;
    # make the content
    my $protocol ||= "acedb";
    my $username ||= "any";
    my $password ||= "any";
    my $host     ||= hostname();
    my $port     ||= $self->{'_SGIF_LOCAL_SERVER'}->port;
    my $seq      ||= 0;
    my $writebck ||= 0;
    my $style    ||= "ZMap.styles";
    my $content    = sprintf($fmt,
                             $protocol,
                             $username,
                             $password,
                             $host,
                             $port,
                             ($seq ? 'true' : 'false'),
                             ($writebck ? 'true' : 'false'),
                             $style);
    return $content;
}



sub make_request{
    my ($self, @commands) = @_;

    my $xr = $self->current_xclient;
    unless($xr){
        $self->write_on_canvas("No current window.");
        return;
    }

    push(@commands, 'zoom_in') unless @commands;
    warn "$self requesting @commands\n";
    my @a = $xr->send_commands(@commands);
    my $canvasMessage = "Sent commands\n";
    for(my $i = 0; $i < @commands; $i++){
        $canvasMessage .= " - cmd: '" . $commands[$i] . "',  result: '" . $a[$i] . "'\n";
        my ($status, $n) = parse_response($a[$i]);
        if($status == 412){
            delete_xclient_with_id($n->{'error'}->{'windowid'});
        }

        if(ref($n->{response}) eq 'HASH'){
            my ($name, $id) = ($n->{response}->{zmapid}, $n->{response}->{windowid});
            if($name){
                xclient_with_name($name, $id) if $id;
                $self->set_entry_value($name);
            }
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
#    if(!$cnnct){
        my $req = $self->zmapConnectObj->connect_request;
        $self->make_request($req);
        $self->{'_connectedTo'}->{$id} = 1;
 #   }
}

#=====================================================================
# RECEIVE_FILTER && functions it calls
#=====================================================================
sub RECEIVE_FILTER{ 
    my ($_connect, $_request, $_obj, @valid_list) = @_; 
    my ($_status, $_response) = (404, "<error>Unknown Command</error>");

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

    my $xr = xclient_with_name('main');
    my $z  = $self->zmapConnectObj;
    my $h  = {
        response => {
            client => [{
                created => 0,
                exists  => 1,
            }]
        }
    };
    $z->protocol_add_meta($h);


    unless($p->{'id'} 
           && $p->{'request'}
           && $p->{'response'}){
        warn "mismatched request for register_client:\n", 
        "id, request and response required\n",
        "Got '$request'\n";
        return (403,$z->basic_error("Bad Request"));
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

    # this feels convoluted
    $h->{'response'}->{'client'}->[0]->{'created'} = 1;

    return (200, make_xml($h));
}
sub watching_sub{
    my ($watch, $value) = @_;
    my ($self) = @{$watch->Args('-store')};
    my $xr     = xclient_with_name('main');
#    my @clones = qw(b0250);
    my @clones = qw(1.141003345-141198932);
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
    $SIG{CHLD} = 'DEFAULT';
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
