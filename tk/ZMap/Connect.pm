package ZMap::Connect;

=pod

=head1 NAME 

ZMap::Connect

=head1 DESCRIPTION

For  connecting  to ZMap  and  makes  the  X11::XRemote module  more
useable as a server.

=cut

use strict;
use warnings;
use X11::XRemote;
use Tk::X;

my $DEBUG_CALLBACK = 0;
my $DEBUG_EVENTS   = 0;

=head1 METHODS

=head2 new([-option => "value"])

Creates a  new ZMap::Connect Object  of the style  requested.  Control
the style with the following options (defaults in []).

=over 10

=item I<-server> '0|1' [1]

 Whether to create a server.

=item I<-request_name> name [_PERL_REQUEST_NAME]

 The atom name.

=item I<-response_name> atom name for response [_PERL_RESPONSE_NAME]

 The atom name.

=item I<-debug> '0|1' [0]

 Show debugging messages when set to 1

=back

=cut

sub new{
    my ($pkg, @args) = @_;
    my $p = {};
    $p = { @args } if(!(scalar(@args) % 2));
    my $self = {_server        => 1,
                _debug         => 0,
                _request_name  => '_PERL_REQUEST_NAME',
                _response_name => '_PERL_RESPONSE_NAME',
                map { "_".lc(substr($_,1)) => $p->{$_} } keys(%$p)
                };
    bless($self, $pkg);
    return $self;
}


=head2 init(Tk, [handler, [data]])

Initialises  the new  object  so  it is  useable.  Requires the  B<Tk>
object.  The  B<handler> is  the callback which  gets called  when the
window is sent a message via an  atom.  It is called by THIS module as
C<<< $callback->($self, $request, @data) >>>.  Note the data supplied as
a  list ref  to this  function gets  dereferenced when  passed  to the
callback.   The request  string  is  supplied free  of  charge so  the
callback does  not have to  navigate its way  to it from  this module
(supplied as $self above).

 Usage:
 my $fruits = [qw(apples pears bananas)];
 my $veg    = [qw(carrots potatoes)];
 my $callback = sub{ 
     my ($zmap, $req, $f, $v) = @_;
     if($req =~ /fruits/){
        print "I know about these fruits: " . join(", ", @$f) . "\n";
     }elsif($req =~ /veg/){
        print "I know about these vegetables: " . join(", ", @$v) . "\n";
     }
     return (200, "printed my knowledge.");
 };
 my $zmap = ZMap::Connect->new();
 $zmap->init($tk, $callback, [$fruits, $veg]);

=cut

sub init{
    my ($self, $tk, $callback, $data) = @_;
    unless($tk){
        warn "usage: ".__PACKAGE__."->init(Tk_Object, [SubRoutine_Ref, Data_Array_Ref]);\n"; 
        return;
    }
    $self->requestWidget($tk);
    $self->responseWidget($tk);
    $self->respond_handler($callback, $data);
}

=head2 connect_request( )

This maybe  used to  get the  string for registering  the server  as a
remote window  of the  zmap client.  ZMap  should honour  this request
creating a  client with  the window  id of the  request widget  as its
remote window.

=cut

sub connect_request{
    my ($self) = @_;
    my $req = $self->requestWidget();
    my $xr  = $self->xremote();
    my $fmt = "%s id = %s ; request = %s ; response = %s ;";
    return sprintf($fmt, 
                   "register_client", 
                   $req->id, 
                   $self->request_name,
                   $self->response_name
                   );
}

=head2 server_window_id( )

Just the window id of the request widget.

=cut


sub server_window_id{
    my ($self) = @_;
    return undef unless $self->_is_server();
    return $self->requestWidget->id();
}

=head2 xremote( )

The xremote Object [C<<< X11::XRemote >>>].

=cut

sub xremote{
    my ($self, $id) = @_;
    my $xr = $self->{'_xremote'};
    if(!$xr){
        if($self->_is_server){
            $xr = X11::XRemote->new(-server => 1,
                                    -id     => $id
                                    );
        }else{
            $xr = X11::XRemote->new();
        }
        $self->{'_xremote'} = $xr;
    }
    return $xr;
}


=head1 SETUP METHODS

It is  recommended that these methods  are not used  directly to setup
the object (as  setters), but feel free to call  with no arguments (as
getters).

=head2 request_name( )

Set/Get the atom name for the request.

=cut

# ======================================================== #
# SETUP STUFF: Easier just to call $self->init(@args);     #
# ======================================================== #
sub request_name{
    my ($self, $name) = @_;
    if($name && !$self->requestWidget()){
        $self->{'_request_name'} = $name;
    }
    return $self->{'_request_name'};
}

=head2 response_name( )

Set/Get the atom name for the response.

=cut

sub response_name{
    my ($self, $name) = @_;
    if($name && !$self->responseWidget()){
        $self->{'_response_name'} = $name;
    }
    return $self->{'_response_name'};
}

=head2 requestWidget( )

Set/Get the request widget.

=cut

sub requestWidget{
    my ($self, $tk) = @_;
    my $label = $self->{'_requestWidget'};
    if($tk && !$label){
        my $aName = $self->request_name();
        $label    = $tk->Label(
                               -text => "${aName}Widget",
                               )->pack(-side => 'left');
        # This actually sets up the XREMOTE.
        my $xr = $self->xremote($label->id());
        #$xr->request_name($label->atomname($label->atom($label->PathName)));
        $xr->request_name($aName);
        $label->packForget();
        $self->{'_requestWidget'} = $label;
    }
    return $label;
}

=head2 responseWidget( )

Set/Get the response widget.

=cut

sub responseWidget{
    my ($self, $tk) = @_;
    my $label = $self->{'_responseWidget'};
    if($tk && !$label){
        my $aName = $self->response_name();
        $label    = $tk->Label(
                               -text => "${aName}Widget"
                               )->pack(-side => 'left');

        my $xr = $self->xremote();
        #$xr->response_name($label->atomname($label->atom($label->PathName)));
        $xr->response_name($aName);

        $label->packForget();
        $self->{'_responseWidget'} = $label;
    }
    return $self->{'_responseWidget'};
}

=head2 respond_handler( )

Set/Get the callback which will get called.

=cut

sub respond_handler{
    my ($self, $callback, $data) = @_;
    $self->__callback($callback);
    $self->__callback_data($data);
    my $handler = \&ZMap::Connect::_do_callback;
    if(my $widget = $self->requestWidget){
        $widget->Tk::bind('<Property>', [ $handler , $self ] );
    }else{
        warn "Sugested usage:\n" . 
            "my \$c = ".__PACKAGE__."->new([options]);\n" .
            "\$c->init(\$tk, \$callback, \$callback_data);\n";
    }
}

# ======================================================== #
#                      INTERNALS                           #
# ======================================================== #
sub _do_callback{
    my ($tk, $self) = @_;
    my $ev    = $tk->XEvent(); # Get the event
    my $state = ($ev->s ? $ev->s : 0); # assume zero (PropertyDelete I think)
    my $reqnm = $self->request_name(); # atom name of the request
    return if $state == PropertyDelete;
    #====================================================================
    # DEBUG STUFF
    warn "\n//========== _doing_callback ==========\n" if $DEBUG_CALLBACK;
    #$state == PropertyNewValue;
    if($DEBUG_EVENTS){
        foreach my $m('a'..'z','A'..'Z'){
            warn "Event on ".($self->is_server ? "server" : "client").
                " method '$m' - ". $ev->$m() . " \n" if $ev->$m();
        }
    }
    unless($ev->T eq 'PropertyNotify'){ warn "Not a propertyNotify\n"; }
    unless($ev->d eq $reqnm){
        warn "Event was NOT for this ".($self->is_server ? "server" : "client").
            "\n" if $DEBUG_CALLBACK;
        return ;
    }
    my $req = $self->xremote->request_string();
    warn "Event has request string $req\n" if $DEBUG_CALLBACK;
    #=========================================================
    my $cb = $self->__callback();
    my @data = @{$self->__callback_data};
    my $reply;
    my $fstr  = $self->xremote->format_string;
    eval{ 
        my ($status,$str) = $cb->($self, $req, @data);
        $reply = sprintf($fstr, $status, sprintf("<xml>%s</xml>", $str));
    };
    if($@){
        $reply ||= sprintf($fstr, 500, "<xml><simplemessage>Internal Server Error $@</simplemessage></xml>");
    }
    $reply ||= sprintf($fstr, 500, "<xml><simplemessage>Internal Server Error</simplemessage></xml>");
    $self->xremote->send_reply($reply);
}

sub __callback_data{
    my($self, $dataRef) = @_;
    $self->{'_callback_data'} = $dataRef if ($dataRef && ref($dataRef) eq 'ARRAY');
    return $self->{'_callback_data'} || [];
}
sub __callback{
    my($self, $codeRef) = @_;
    $self->{'_callback'} = $codeRef if ($codeRef && ref($codeRef) eq 'CODE');
    return $self->{'_callback'} || sub { warn "@_\nNo callback set.\n"; return (500,"") };
}
sub _is_server{
    my ($self) = @_;
    return ($self->{'_server'} ? 1 : 0);
}
# ======================================================== #
# DESTROY: Hopefully won't need to do anything             #
# ======================================================== #
sub DESTROY{
    my ($self) = @_;
    warn "Destroying";
}

1;

__END__

=head1 AUTHOR

R.Storey <rds@sanger.ac.uk>

=head1 SEE ALSO

L<X11::XRemote>, L<perl>.

=cut
