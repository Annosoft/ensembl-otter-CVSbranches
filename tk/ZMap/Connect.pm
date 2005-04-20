package ZMap::Connect;

=pod

=head1 NAME 

 ZMap::Connect

=head1 DESCRIPTION

 For  connecting  to ZMap  and  makes  the  X11::XRemote modules  more
 useable as a server.

=cut

use strict;
use warnings;
use X11::XRemote;
use Tk::X;

=head1 METHODS

=head2 new()

 Creates a new ZMap::Connect Object.

=cut

my $DEBUG_CALLBACK = 0;
my $DEBUG_EVENTS   = 0;

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


=head1 init(Tk, handler, data)

 Initialises the new object so it is useable. Requires the Tk object.

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

sub server_window_id{
    my ($self) = @_;
    return undef unless $self->_is_server();
    return $self->requestWidget->id();
}

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

sub _is_server{
    my ($self) = @_;
    return ($self->{'_server'} ? 1 : 0);
}

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
sub response_name{
    my ($self, $name) = @_;
    if($name && !$self->responseWidget()){
        $self->{'_response_name'} = $name;
    }
    return $self->{'_response_name'};
}
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

# ======================================================== #
# DESTROY: Hopefully won't need to do anything             #
# ======================================================== #
sub DESTROY{
    my ($self) = @_;
    warn "Destroying";
}

1;
