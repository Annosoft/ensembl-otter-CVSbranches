package CanvasWindow::XRemoteExample3;

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
use Data::Dumper;

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
        reset_sigCHLD(); # don't need to catch anymore
        $self = undef;
    };
    $self->canvas->toplevel->protocol('WM_DELETE_WINDOW',  $close_window);
    # Create a Button Frame for the buttons, etc...
    my $button_frame2 = $self->canvas->toplevel->Frame->pack(-side => 'top', -fill => 'x');
    my $button_frame1 = $self->canvas->toplevel->Frame->pack(-side => 'top', -fill => 'x');
    # Add an Entry widget for the window ID


    my $kill = $button_frame1->Button(-text => 'End Test',
                                     -command => sub {
                                         my $top = $self->canvas->toplevel;
                                         $top->destroy;
                                     }
                                     )->pack(-side => 'left');


    my $idEntry = $button_frame2->BrowseEntry(-label    => "Window ID:",
                                             -width    => 20,
                                             -variable => $self->win_id_entry_ref(),
                                             -font     => ['Helvetica', 12, 'normal'],
                                             -listcmd  => sub {
                                                 my $w = shift;
                                                 $w->choices([ ZMap::ConnectUtils::list_xclient_ids() ])
                                                 },
                                             )->pack(-side => 'left');

    my $query = sprintf('feature_find type = %s ; method = %s ; start = %d ; end = %d ; q_start = %d ; q_end = %d ; strand = %s ; feature = %s',
                        qw(homol wublastx_worm 10788 11048 514 600 + WP:CE19655));
    
    $self->set_entry_value('req_entry_ref', $query);
    my $requestEntry = $button_frame2->Entry(
                                             -width    => 80,
                                             -textvariable => $self->req_entry_ref(),
                                             -font     => ['Helvetica', 12, 'normal'],
                                             )->pack(-side => 'left');
    # Setup some buttons to do things.
    my $open = 
        $button_frame2->Button(-text       => 'Do',
                              -command    => sub { 
                                  $self->make_Entry_request();
                              })->pack(-side => 'right');
    my $zoomIn = 
        $button_frame1->Button(-text       => 'Zoom In',
                              -command    => sub { 
                                  $self->make_request('zoom_in');
                              })->pack(-side => 'right');
    my $zoomOut = 
        $button_frame1->Button(-text       => 'Zoom Out',
                              -command    => sub { 
                                  $self->make_request('zoom_out');
                              })->pack(-side => 'right');


    return $self;
}

#=====================================================================
# Button Functions
#=====================================================================

sub make_Entry_request{
    my ($self) = @_;
    my $request = $self->req_entry_ref();
    $self->make_request($$request);
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
        $canvasMessage .= " - cmd: '" . $commands[$i] . "',\n - result: '" . $a[$i] . "'\n";
        my ($status, $n) = parse_response($a[$i]);
        my ($name, $id) = ($n->{zmapid}, $n->{windowid});
    }
    $self->write_on_canvas($canvasMessage);
}

#=====================================================================
# short helper bits
#=====================================================================
sub current_xclient{
    my ($self) = @_;
    return xclient_with_id(${$self->win_id_entry_ref});
}
sub win_id_entry_ref{
    my ($self) = @_;
    my $n = '';
    my $k = (caller(0))[3];
    $self->{$k} ||= \$n;
    return $self->{$k};
}
sub req_entry_ref{
    my ($self) = @_;
    my $n = '';
    my $k = (caller(0))[3];
    $self->{$k} ||= \$n;
    return $self->{$k};
}
sub set_entry_value{
    my ($self, $name, $value) = @_;
    my $ref = $self->$name();
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
