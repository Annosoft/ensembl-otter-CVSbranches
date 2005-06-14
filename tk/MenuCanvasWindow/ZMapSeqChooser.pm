### MenuCanvasWindow::ZmapSeqChooser

package MenuCanvasWindow::ZMapSeqChooser;

use strict;
use Carp qw{ cluck confess };
use MenuCanvasWindow::SeqChooser;
use ZMap::Connect qw(:all);
use Sys::Hostname;
use Tie::Watch;

our @ISA = qw(MenuCanvasWindow::SeqChooser);

sub update_display{
    my ($self , $ace) = @_ ;
    warn "cannot update display yet";
    return 1;
}
# this should be called when a user tries to save, but no ZMap is opened
sub open_dialogue{
    my ($self) = @_ ;
    warn "cannot open dialogue yet";
    return undef;

}
sub attach {
    my( $self ) = @_;
    
    if (my $xwid = $self->zmap_id) {
        my $xrem = Hum::Ace::XaceRemote->new($xwid);
        $self->xace_remote($xrem);
        #$xrem->send_command('save');
# This command may be redundant with zmap + lace
#        $xrem->send_command('writeaccess -gain');
    } else {
        warn "no xwindow id: $xwid";
    }
}

sub process_id_list {
    my( $self, $zmap_process_id ) = @_;
    
    if ($zmap_process_id) {
        $self->{'_zmap_process_id'} = $zmap_process_id;
    }
    return $self->{'_zmap_process_id'};
}
sub launch {
    my( $self ) = @_;
    
    $self->kill;
    my $z = $self->insert_zmap_connector();
    my $s_id = $z->server_window_id();
    if (my $path = $self->ace_path) {
        if($self->spawn_sgifaceserver){
            $self->write_dot_zmap();
            sleep(2);
            if (my $pid = fork) {
                print STDERR " zmap forked\n";
                $self->process_id_list($pid);
            }
            elsif (defined($pid)) {
                my $c = $self->zmap_dir();
                my @e = ('zmap', '--conf_dir', $c, '--win_id', $s_id);
                print  "I ran '@e' \n";
                exec( @e );
            }
            else {
                confess "Error: can't fork : $!";
            }

        } else {
             warn "Error: couldn't spawn sgifaceserver for some reason...";
        }
    } else {
        warn "Error: ace_path not set";
    }
}
sub kill {
    my( $self ) = @_;
    
    if (my $pid = $self->process_id_list) {
        warn "Killing zmap process '$pid'\n";
        kill 9, $pid;
    }
}

sub remote {}

sub get_window_id {
    my( $self ) = @_;

    # be good if we could replace this with something more automatic....    
    my $mid = $self->message("Please click on the xace main window with the cross-hairs");
    $self->delete_message($mid);
    local *XWID;
    open XWID, "xwininfo |"
        or confess("Can't open pipe from xwininfo : $!");
    my( $xwid );
    while (<XWID>) {
        # xwininfo: Window id: 0x7c00026 "ACEDB 4_9c, lace bA314N13"

      # HACK
      # above format NOT returnd by xwininfo on Sun OS 5.7:
      #  tace version:
      #  ACEDB 4_9r,  build dir: RELEASE.2003_05_01
      # 2 lines before modified to support xace at RIKEN

        # BEFORE: if (/Window id: (\w+) "([^"]+)/) {
        if (/Window id: (\w+) "([^"]+)/ || /Window id: (\w+)/) {
            my $id   = $1;
            my $name = $2;
	    # BEFORE: if ($name =~ /^ACEDB/){
            if ($name =~ /^ACEDB/ || $name eq '') {
                $xwid = $id;
                $self->message("Attached to:\n$name");
            } else {
                $self->message("'$name' is not an xace main window");
            }
        }
    }
    if (close XWID) {
        return $xwid;
    } else {
        $self->message("Error running xwininfo: $?");
    }
}


sub DESTROY {
    my( $self ) = @_;
    # need to undef AceDatabase for it's clean up.
    # warn "AceDatabase->unlock should now happen\n";
    $self->kill(); # get rid of zmap
    delete $self->{'_SGIF_LOCAL_SERVER'}; # shutdown server
    delete $self->{'_AceDatabase'}; # unlock will now happen
    # warn "AceDatabase->unlock should have happened\n";
    if (my $sn = $self->SequenceNotes){
        # then refresh locks
        # warn "lock refresh should now happen\n";
        $sn->refresh_column(7) ; ## locks column
        # warn "lock refresh should have happened\n";
        # need to clean up the sequenceNotes reference
        delete $self->{'_sequence_notes'} ;
    }
    warn "Destroying ZmapSeqChooser for ", $self->ace_path, "\n";
}

#==============================================================================#
#
# EXTRA METHODS...
#
#==============================================================================#

sub insert_zmap_connector{
    my ($self) = @_;
    my $justOneInstance = 
    my $zc = $self->{'_zmap_connector'};
    if(!$zc){
        my $zmap = ZMap::Connect->new( -server => 1 );
        $zmap->init($self->menu_bar(), \&RECEIVE_FILTER, [ $self, qw( register_client ) ]);
        my $id = $zmap->server_window_id();
        $zc = $self->{'_zmap_connector'} = $zmap;
    }
    return $zc;
}

sub spawn_sgifaceserver{
    my ($self) = @_;
    my $sgif   = $self->{'_SGIF_LOCAL_SERVER'};
    my $path   = $self->ace_path();
    unless($sgif){
        $sgif = Hum::Ace::LocalServer->new($path);
        $sgif->server_executable("/nfs/team71/acedb/edgrif/acedb/CODE/acedb/bin.ALPHA_5/sgifaceserver");
        $sgif->timeout_string('0:0');
        $sgif->start_server() or return 0; # this only check the fork was successful
        $sgif->ace_handle     or return 0; # this checks it can connect
        $self->zmap_port($sgif->port);
        $self->{'_SGIF_LOCAL_SERVER'} = $sgif;
    }
    return 1;
}

sub connector{
    my ($self, $c) = @_;
    $self->{'_zmap_connect'} = $c if $c;
    return $self->{'_zmap_connect'};
}

sub write_dot_zmap{
    my ($self) = @_;
    my $dir    = $self->zmap_dir();
    my $file   = "${dir}/ZMap";
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
    my $fmt = qq`source\n{\nurl = "%s://%s:%s@%s:%d"\nsequence = %s\nwriteback = %s\nstylesfile = "%s"\n}\n`;
    # make the content
    my $protocol ||= "acedb";
    my $username ||= "any";
    my $password ||= "any";
    my $host     ||= hostname();
    my $port     ||= $self->zmap_port;
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

sub zmap_port{
    my ($self, $port)  = @_;
    $self->{'_zmap_port'} = $port if $port;
    return $self->{'_zmap_port'};
}

sub zmap_dir {
    my( $self, $zmap_dir ) = @_;
    warn "Set using the Config file please.\n" if $zmap_dir;
    my $ace_path = $self->ace_path();
    return $self->Client->option_from_array([qw( zmap zmap_dir )]) || "$ace_path/.ZMap";
}



sub RECEIVE_FILTER{
    my ($_connect, $_request, $_obj, @list) = @_;
    my ($_status, $_response) 
        = (404, "<error>unknown command</error>");

    foreach my $valid(@list){
        if($_request =~ s/^$valid//){
            ($_status, $_response) 
                = $_obj->$valid($_request) if $_obj->can($valid);
            last;
        }
    }

    return ($_status, $_response) ;    
}

sub register_client{
    my ($self, $request) = @_;

    my $p  = parse_params($request);
    my $xr = xclient_with_name('main');

    return if $xr;

    unless($p->{'id'} 
           && $p->{'request'}
           && $p->{'response'}){
        warn "mismatched request for register_client:\n", 
        "id, request and response required\n",
        "Got '$request'\n";
        return (403, "<error>Bad request!</error>");
    }
 
    $xr = xclient_with_name('main', $p->{'id'});

    $self->set_entry_value('main');

    Tie::Watch->new(-variable => \$WAIT_VARIABLE,
                    -debug    => 1,
                    -store    => [ \&open_clones, $self ],
                    );

    return (200, "<error>Bingo.</error>");
}

sub open_clones{
    my ($watch) = @_;
    my ($self)  = @{$watch->Args('-store')};

    my @clones = $self->clone_list;
    my ($chr, $st, $end) = split(/\.|\-/, $clones[0]);
    my $cmd = "newZmap seq = @clones ; start = 1 ; end = 0";

    $self->make_request($cmd);
    
    $watch->Unwatch;
}

sub make_request{
    my ($self, @commands) = @_;

    my $xr = $self->current_xclient;
    unless($xr){
        warn "No current window.";
        return;
    }

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
    warn $canvasMessage;
}

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


1;

__END__

=head1 NAME - MenuCanvasWindow::ZmapSeqChooser

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

__DATA__
