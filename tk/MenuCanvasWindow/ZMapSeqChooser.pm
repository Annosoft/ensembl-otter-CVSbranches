### MenuCanvasWindow::ZmapSeqChooser

package MenuCanvasWindow::ZMapSeqChooser;

use strict;
use Carp qw{ cluck confess };
use MenuCanvasWindow::SeqChooser;
use ZMap::Connect;
use Sys::Hostname;
use Tie::Watch;

our @ISA = qw(MenuCanvasWindow::SeqChooser);

sub insert_zmap_connector{
    my ($self) = @_;
    my $justOneInstance = '_zmap_connector';
    my $zc = $self->{$justOneInstance};
    if(!$zc){
        my $zmap = ZMap::Connect->new( -server => 1 );
        $zmap->init($self->menu_bar(), \&RECEIVE_FILTER, [ $self, qw( register_client ) ]);
        my $id = $zmap->server_window_id();
        $self->{'_xclients'}->{$id} = $zmap->xremote;
        $zc = $self->{$justOneInstance} = $zmap;
    }
    return $zc;
}

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
sub xclient_with_name{
    my ($self, $name, $id) = @_;
    if(!$id){
        $id = $self->{$name};
    }else{ 
        $self->{$name} = $id;
    }
    return unless $id;
    return $self->xclient($id);
}
sub update_display{
    my ($self , $ace) = @_ ;
    my $xr = $self->remote  || $self->open_dialogue;
    print STDERR "Sending:\n$ace";
    if ($xr) {
         # I think this lot can be replaced with a single command to zmap.
#        $xr->load_ace($ace);
#        $xr->save;
#        $xr->send_command('gif ; seqrecalc');
        $xr->send_command('edit  ***should send details of thing to be recalculated here');
        return 1;
    } else {
        $self->message("No zmap attached");
        print STDERR "Not able to update sequence display - no zmap attached";
        return 0;
    }    
}

# this should be called when a user tries to save, but no ZMap is opened
sub open_dialogue{
    my ($self) = @_ ;
    
    my $answer = $self->canvas->toplevel()->messageBox(-title => 'Please Reply', 
     -message => 'No ZMap attached, would you like to launch ZMap ?', 
     -type => 'YesNo', -icon => 'question', -default => 'Yes');

    if ($answer eq 'Yes'){
        $self->launch_zmap();
    }
    return $self->xace_remote;
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

# This will all need changing.........
sub get_xwindow_id_from_readlock__________________________________________________________ {
    my( $self ) = @_;
    
    local(*ZMAP_DIR, *ZMAP_FILE);

    my $pid  = $self->process_id_list or confess "zmap_process_id not set";
    my $path = $self->zmap_dir        or confess "zmap_dir not set";
    
    # Find a zmap windowid file in the zmap dir for the process we just launched
    my $zmap_dir = $path;
    my( $zmap_file );
    my $wait_seconds = 20;
    for (my $i = 0; $i < $wait_seconds; $i++, sleep 1) {
        opendir ZMAP_DIR, $zmap_dir or confess "Can't opendir '$zmap_dir' : $!";

        print STDERR "I'm looking in $zmap_dir";

        ($zmap_file) = grep /main\.win_id/, readdir ZMAP_DIR;

        print STDERR "hmmm $zmap_file";

        closedir ZMAP_DIR;

        if ($zmap_file) {
            $zmap_file = "$zmap_dir/$zmap_file";
            last;
        }
    }
    unless ($zmap_file) {
        warn "Can't find zmap window id file in '$zmap_dir' for process '$pid' after waiting for $wait_seconds seconds\n";
        return 0;
    }
    
    my( $xwid );
    for (my $i = 0; $i < $wait_seconds; $i++, sleep 1) {
        # Extract the WindowID from the zmap file
        open ZMAP_FILE, $zmap_file or confess "Can't read '$zmap_file' : $!";
        while (<ZMAP_FILE>) {
            #warn "Looking at: $_";
            if (/WindowID: (\w+)/) {
                $xwid = $1;
                last;
            }
        }
        close ZMAP_FILE;
        
        last if $xwid;
    }
    if ($xwid) {
        $self->xclient_with_name('main', $xwid);

        return 1;
    } else {
        warn "WindowID was not found in zmap file - did zmap start up correctly ?";
        return 0;
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

#==============================================================================#
#
# EXTRA METHODS...
#
#==============================================================================#


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




sub register_client{
    my ($self, $request) = @_;

    my $p  = parse_params($request);
    my $xr = $self->xclient_with_name('main');

    return if $xr;

    unless($p->{'id'} 
           && $p->{'request'}
           && $p->{'response'}){
        warn "mismatched request for register_client:\n", 
        "id, request and response required\n",
        "Got '$request'\n";
        return (403, "<error>Bad request!</error>");
    }
 
    $xr = $self->xclient_with_name('main', $p->{'id'});
    print "Here";

    Tie::Watch->new(-variable => \$WAIT_VARIABLE,
                    -debug    => 1,
                    -store    => [ \&open_clones, $self ],
                    );

    return (200, "<error>Bingo.</error>");
}

sub open_clones{
    my ($watch) = @_;
    my ($self)  = @{$watch->Args('-store')};
    my $xr = $self->xclient_with_name('main');

    my @clones = $self->clone_list;
    my ($chr, $st, $end) = split(/\.|\-/, $clones[0]);
    my $cmd = "newZmap seq = @clones ; start = 1 ; end = 0";
    print "I'm sending $cmd\n";
    $xr->send_commands($cmd);
    
    $watch->Unwatch;
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

1;

__END__

=head1 NAME - MenuCanvasWindow::ZmapSeqChooser

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

__DATA__
