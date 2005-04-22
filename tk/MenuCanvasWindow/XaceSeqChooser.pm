
### MenuCanvasWindow::XaceSeqChooser

package MenuCanvasWindow::XaceSeqChooser;

use strict;
use warnings;
use Hum::Ace::XaceRemote;
use Carp qw{ cluck confess };

use MenuCanvasWindow::SeqChooser;
our @ISA = qw(MenuCanvasWindow::SeqChooser);

sub update_display{
    my ($self , $ace) = @_ ;

    my $xr = $self->remote  || $self->open_dialogue;
   
    print STDERR "Sending:\n$ace";
    if ($xr) {
        $xr->load_ace($ace);
        $xr->save;
        $xr->send_command('gif ; seqrecalc');
       
        return 1;
    } else {
        $self->message("No xace attached");
        print STDERR "not able to send .ace file - no xace attached";
        return 0;
    }    
}

# this should be called when a user tries to save, but no Xace is opened
sub open_dialogue{
    my ($self) = @_ ;
    
    my $answer = $self->canvas->toplevel()->messageBox(-title => 'Please Reply', 
     -message => 'No Xace attached, would you like to launch Xace?', 
     -type => 'YesNo', -icon => 'question', -default => 'Yes');

    if ($answer eq 'Yes'){
        $self->launch();
    }
    return $self->remote;
}

sub attach {
    my( $self ) = @_;
    
    if (my $xwid = $self->get_window_id) {
        my $xrem = Hum::Ace::XaceRemote->new($xwid);
        $self->remote($xrem);
        #$xrem->send_command('save');
        $xrem->send_command('writeaccess -gain');
    } else {
        warn "no xwindow id: $xwid";
    }
}

sub process_id_list {
    my( $self, $xace_process_id ) = @_;
    
    if ($xace_process_id) {
        $self->{'_xace_process_id'} = $xace_process_id;
    }
    return $self->{'_xace_process_id'};
}

sub launch {
    my( $self ) = @_;
    
    $self->kill;
    
    if (my $path = $self->ace_path) {
        if (my $pid = fork) {
            $self->process_id_list($pid);
            $self->get_xwindow_id_from_readlock;
        }
        elsif (defined($pid)) {
            exec('xace', '-fmapcutcoords', $path);
        }
        else {
            confess "Error: can't fork : $!";
        }
    } else {
        warn "Error: ace_path not set";
    }
}

sub kill {
    my( $self ) = @_;
    
    if (my $pid = $self->process_id_list) {
        warn "Killing xace process '$pid'\n";
        kill 9, $pid;
    }
}

sub get_xwindow_id_from_readlock {
    my( $self ) = @_;
    
    local(*LOCK_DIR, *LOCK_FILE);

    my $pid  = $self->process_id_list or confess "xace_process_id not set";
    my $path = $self->ace_path        or confess "ace_path not set";
    
    # Find the readlock file for the process we just launched
    my $lock_dir = "$path/database/readlocks";
    my( $lock_file );
    my $wait_seconds = 20;
    for (my $i = 0; $i < $wait_seconds; $i++, sleep 1) {
        opendir LOCK_DIR, $lock_dir or confess "Can't opendir '$lock_dir' : $!";
        ($lock_file) = grep /\.$pid$/, readdir LOCK_DIR;
        closedir LOCK_DIR;
        if ($lock_file) {
            $lock_file = "$lock_dir/$lock_file";
            last;
        }
    }
    unless ($lock_file) {
        warn "Can't find xace readlock file in '$lock_dir' for process '$pid' after waiting for $wait_seconds seconds\n";
        return 0;
    }
    
    my( $xwid );
    for (my $i = 0; $i < $wait_seconds; $i++, sleep 1) {
        # Extract the WindowID from the readlock file
        open LOCK_FILE, $lock_file or confess "Can't read '$lock_file' : $!";
        while (<LOCK_FILE>) {
            #warn "Looking at: $_";
            if (/WindowID: (\w+)/) {
                $xwid = $1;
                last;
            }
        }
        close LOCK_FILE;
        
        last if $xwid;
    }
    if ($xwid) {
        my $xrem = Hum::Ace::XaceRemote->new($xwid);
        $self->remote($xrem);
        #$xrem->send_command('save');
        $xrem->send_command('writeaccess -gain');
        return 1;
    } else {
        warn "WindowID was not found in lock file - outdated version of xace?";
        return 0;
    }
}

sub remote {
    my( $self, $xrem ) = @_;
    
    if ($xrem) {
        my $expected = 'Hum::Ace::XaceRemote';
        confess "'$xrem' is not an '$expected'"
            unless (ref($xrem) and $xrem->isa($expected));
        $self->{'_xace_remote'} = $xrem;
    }
    return $self->{'_xace_remote'};
}

sub get_window_id {
    my( $self ) = @_;
    
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
    warn "Destroying XaceSeqChooser for ", $self->ace_path, "\n";
}
