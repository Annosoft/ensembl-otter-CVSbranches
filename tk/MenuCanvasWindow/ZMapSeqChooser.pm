### MenuCanvasWindow::ZmapSeqChooser

package MenuCanvasWindow::ZMapSeqChooser;

use strict;
use Carp qw{ cluck confess };
use MenuCanvasWindow::SeqChooser;
use ZMap::Connect qw(:all);
use Sys::Hostname;
use Tie::Watch;
use Hum::Ace::LocalServer;

my $ZMAP_DEBUG = 1;
#==============================================================================#
#
# WARNING: THESE ARE INJECTED METHODS!!!!
#  I HAVE PREFIXED THEM ALL WITH zMap SO NONE SHOULD CLASH
#  BUT ALL WILL NEED CHANGING LATER (RDS)
#
#==============================================================================#

=pod

=head1 WARNING

This modules injects methods into the MenuCanvasWindow::XaceSeqChooser
namespace.  All have  been prefixed with "zMap" to  avoid any clashes,
but this isn't a long term solution.

=head1 zMapLaunchZmap

This is where it all starts.  This is the method which gets called 
on 'Launch ZMap' menu item in xaceseqchooser window.

=cut

sub MenuCanvasWindow::XaceSeqChooser::zMapLaunchZmap {
    my( $self ) = @_;

    $self->zMapKillZmap;
    if (my $path = $self->ace_path) {
        if(my $ok = $self->zMapSpawnSgifaceserver){
            my $z = $self->zMapInsertZmapConnector();
            $self->zMapWriteDotZmap();
            $self->isZMap(1);
            my @e = ('zmap', 
                     '--conf_dir' => $self->zMapZmapDir,
                     '--win_id'   => $z->server_window_id);
            warn "export PATH=$ENV{PATH}\nexport LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}\n@e\n" if $ZMAP_DEBUG;
            sleep(2);
            my $ref = Hum::Ace::LocalServer::full_child_info();
            my $pid = fork_exec(\@e, $ref, 0, sub { 
                my ($info) = @_;
                flush_bad_windows;
                #use Data::Dumper;                
                #warn "INFO: ", Dumper($info);
            });

            if($pid){
                $self->zMapProcessIDList($pid);
            }else{
                my $mess = "Error: couldn't fork()\n";
                warn $mess;
                $self->message($mess);
            }
        } else {
            my $mess = "Error: couldn't spawn sgifaceserver for some reason...";
            warn $mess;
            $self->message($mess);
        }
    } else {
        warn "Error: ace_path not set";
    }
}
sub MenuCanvasWindow::XaceSeqChooser::zMapKillZmap {
    my( $self ) = @_;
    
    if (my $pid = $self->zMapProcessIDList) {
        warn "Killing zmap process '$pid'\n";
        CORE::kill 9, $pid;
    }
}
sub MenuCanvasWindow::XaceSeqChooser::zMapProcessIDList {
    my( $self, $zmap_process_id ) = @_;
    
    if ($zmap_process_id) {
        $self->{'_zMap_ZMAP_PROCESS_ID'} = $zmap_process_id;
    }
    return $self->{'_zMap_ZMAP_PROCESS_ID'};
}
sub MenuCanvasWindow::XaceSeqChooser::zMapInsertZmapConnector{
    my ($self) = @_;
    my $zc = $self->{'_zMap_ZMAP_CONNECTOR'};
    if(!$zc){
        my $mb   = $self->menu_bar();
        my $zmap = ZMap::Connect->new( -server => 1 );
        $zmap->init($mb, \&RECEIVE_FILTER, [ $self, qw( register_client ) ]);
        my $id = $zmap->server_window_id();
        $zc = $self->{'_zMap_ZMAP_CONNECTOR'} = $zmap;
    }
    return $zc;
}
sub MenuCanvasWindow::XaceSeqChooser::zMapZmapConnector{
    return shift->zMapInsertZmapConnector(@_);
}
sub MenuCanvasWindow::XaceSeqChooser::zMapSpawnSgifaceserver{
    my ($self) = @_;
    my $sgif   = $self->{'_zMap_SGIF_LOCAL_SERVER'};
    my $path   = $self->ace_path();
    unless($sgif){
        $sgif = Hum::Ace::LocalServer->new($path);
        $sgif->server_executable('sgifaceserver');
        $sgif->timeout_string('0:0');
        $sgif->start_server() or return 0; # this only check the fork was successful
        $sgif->ace_handle(1)  or return 0; # this checks it can connect
        $self->zMapPort($sgif->port);
        $self->{'_zMap_SGIF_LOCAL_SERVER'} = $sgif;
    }
    return 1;
}

sub MenuCanvasWindow::XaceSeqChooser::zMapWriteDotZmap{
    my ($self) = @_;
    my $dir    = $self->zMapZmapDir();
    my $file   = "${dir}/ZMap";
    my $fh;
    eval{
        # directory should be made already
        open($fh, ">$file") or die "write_dot_zmap: error writing file '$file', $!";
    };
    warn "Error in :$@" if $@;
    unless($@){
        my $content = $self->zMapDotZmapContent();
        print $fh $content;
        return 1;
    }
    close $fh;
    return 0;
}

sub MenuCanvasWindow::XaceSeqChooser::zMapDotZmapContent{
    my ($self) = @_;
#    my $fmt = qq`source\n{\nurl = "%s://%s:%s@%s:%d"\nsequence = %s\nwriteback = %s\nstylesfile = "%s"\n}\n`;
#    my $fmt = qq`source\n{\nurl = "%s://%s:%s@%s:%d"\nsequence = %s\nwriteback = %s\n}\n`;
    my $fmt = qq`source\n{\nurl = "%s://%s:%s@%s:%d"\nsequence = %s\nwriteback = %s\nfeaturesets = "%s"\n}\n`;
    # make the content
    my $protocol ||= "acedb";
    my $username ||= "any";
    my $password ||= "any";
    my $host     ||= hostname();
    my $port     ||= $self->zMapPort;
    my $seq      ||= 0;
    my $writebck ||= 0;
    my $style    ||= "";#"ZMap.styles";
    # for now we'll get the list of method names and write them into featuresets.
    # the featuresets list controls the order of columns/style and so we get them
    # ordered by right priority...
    my @methods    = $self->zMapListMethodNames_ordered();
    my $sets     ||= join(" ", @methods);
    if(0){
        my @l = qw(Transcript Putative Processed_pseudogene Unprocessed_pseudogene
                   RepeatMasker RepeatMasker_LINE RepeatMasker_SINE
                   trf Predicted_CpG_island 
                   vertebrate_mRNA EST_Human EST_Mouse EST_Other EST
                   BLASTX ensembl Fgenesh Genscan
                   genomewise REFSEQ 
                   Saturated_BLASTX Saturated_EST
                   Saturated_EST_Human Saturated_vertebrate_mRNA);
    }
    my $content    = sprintf($fmt,
                             $protocol,
                             $username,
                             $password,
                             $host,
                             $port,
                             ($seq ? 'true' : 'false'),
                             ($writebck ? 'true' : 'false'),
                             $style || $sets);
    return $content;
}

sub MenuCanvasWindow::XaceSeqChooser::zMapPort{
    my ($self, $port)  = @_;
    $self->{'_zmap_port'} = $port if $port;
    return $self->{'_zmap_port'};
}

sub MenuCanvasWindow::XaceSeqChooser::zMapZmapDir {
    my( $self, $zmap_dir ) = @_;
    warn "Set using the Config file please.\n" if $zmap_dir;
    my $ace_path = $self->ace_path();
    my $path = "$ace_path/.ZMap";
    unless(-d $path){
        mkdir($path, 0777);
        die "Can't mkdir('$path') : $!\n" unless -d $path;
    }
    return $path;
}

sub MenuCanvasWindow::XaceSeqChooser::zMapListMethodNames_ordered{
    my $self = shift;
    my @list = ();
    my $collection = $self->AceDatabase->get_default_MethodCollection;
    foreach my $method (sort {$a->right_priority <=> $b->right_priority}
                        @{$collection->get_all_Methods}) {
        push(@list, $method->name);
    }
    return @list
}

#===========================================================

sub MenuCanvasWindow::XaceSeqChooser::zMapCurrentXclient{
    my ($self) = @_;
    return xclient_with_name(${$self->zMapEntryRef});
}
sub MenuCanvasWindow::XaceSeqChooser::zMapEntryRef{
    my ($self) = @_;
    my $n = '';
    $self->{'_zMap_ENTRY_REF'} ||= \$n;
    return $self->{'_zMap_ENTRY_REF'};
}
sub MenuCanvasWindow::XaceSeqChooser::zMapSetEntryValue{
    my ($self, $value) = @_;
    my $ref = $self->zMapEntryRef();
    $$ref   = $value;
}


sub MenuCanvasWindow::XaceSeqChooser::zMapRegisterClient{
    my ($self, $request) = @_;

    my $p  = parse_request($request);
    my $xr = xclient_with_name('main');
    my $z  = $self->zMapZmapConnector();
    my $h  = {
        response => {
            client => [{
                created => 0,
                exists  => 1,
            }]
        }
    };
    $z->protocol_add_meta($h);

    $self->zMapEntryRef();
    
    return (200, make_xml($h)) if $xr;

    unless($p->{'client'}->{'xwid'} 
           && $p->{'client'}->{'request_atom'}
           && $p->{'client'}->{'response_atom'}){
        warn "mismatched request for register_client:\n", 
        "id, request and response required\n",
        "Got '$request'\n";
        return (403, $z->basic_error("Bad Request!"));
    }
 
    xclient_with_name('main', $p->{'client'}->{'xwid'});

    $self->zMapSetEntryValue('main');

    Tie::Watch->new(-variable => \$WAIT_VARIABLE,
                    -debug    => 1,
                    -store    => [ \&open_clones, $self ],
                    );
    # this feels convoluted
    $h->{'response'}->{'client'}->[0]->{'created'} = 1;
    return (200, make_xml($h));
}


#===========================================================

sub MenuCanvasWindow::XaceSeqChooser::zMapMakeRequest{
    my ($self, $xmlObject, $action) = @_;

    my $xr = $self->zMapCurrentXclient;
    unless($xr){
        warn "No current window.";
        return ;
    }
    warn "Current window " . $xr->window_id . " @_\n" if $ZMAP_DEBUG;
    my $handler ||= \&RESPONSE_HANDLER;
    my $error   ||= \&ERROR_HANDLER;

    my @commands = obj_make_xml($xmlObject, $action);
    warn "@commands" if $ZMAP_DEBUG;
    my @a = $xr->send_commands(@commands);

    for(my $i = 0; $i < @commands; $i++){
        my ($status, $xmlHash) = parse_response($a[$i]);
        if($status =~ /^2\d\d/){ # 200s
            $handler->($self, $xmlHash);
        }else{
            $error->($self, $status, $xmlHash);
        }
    }
    return ;
}

sub MenuCanvasWindow::XaceSeqChooser::zMapUpdateMenu{
    my ($self) = @_;
    
    my $menu_item = $self->{'_zMapMenuItem'};

    my $cleanUpMenu = sub{
        my ($menuRoot, $this) = @_;
        my @current = list_xclient_names();
        $this->{'_zMapSubMenuItems'} ||= {};
        my @remove = ();
        foreach my $k(keys(%{$this->{'_zMapSubMenuItems'}})){
            my $idx = $this->{'_zMapSubMenuItems'}->{$k};
            print "This is $k \n" if $ZMAP_DEBUG;
            if(!(grep /^$k$/, @current)){
                delete $this->{'_zMapSubMenuItems'}->{$k};
                push(@remove, $idx);
            }
        }
        map { $menuRoot->delete($_) } @remove;
    };

    my $fullCleanUpMenu = sub {
        my ($button, $this) = @_;
        if(my $menu = $button->cget('-menu')){
            $cleanUpMenu->($menu, $this);
        }        
    };

    unless($menu_item){
        my $frame  = $self->menu_bar;
        my $button = $frame->Menubutton(-text => 'ZMap')->pack(-side => 'left');
        my $menu   = $button->Menu(-tearoff => 0);
        $button->configure(-menu => $menu);
        $button->bind('<Button-1>', [ $fullCleanUpMenu, $self ]);

        $self->{'_zMapMenuItem'} = 
            $menu_item = $button; #$self->make_menu('ZMap', 1);
    }

    my $addSubMenuItem = sub {
        my ($this, $menuRoot, $name) = @_;
        $this->{'_zMapSubMenuItems'} ||= {};
        $this->{'_zMapSubMenuItems'}->{$name} = scalar(keys(%{$this->{'_zMapSubMenuItems'}}));
        my $submi = $menuRoot->cascade(-label => $name, -tearoff => 0);

        if($name eq 'main'){
            $submi->command(-command => [ sub { 
                my $self = shift;
                my $z = newXMLObj('client');
                $self->zMapSetEntryValue($name); 
                $self->zMapMakeRequest($z, 'finish');
            }, $self ],-label => 'Shutdown ZMap');     
        }else{
            $submi->command(-command => [ sub { 
                my $this = shift;
                my $z    = newXMLObj('zoom_in');
                $this->zMapSetEntryValue($name); 
                $this->zMapMakeRequest($z, 'zoom_in');
            }, $self ], 
                            -label => 'Zoom In');             
            $submi->command(-command => [ sub { 
                my $this = shift;
                my $z    = newXMLObj('zoom_in');
                $self->zMapSetEntryValue($name); 
                $self->zMapMakeRequest($z, 'zoom_out');
            }, $self ], 
                            -label => 'Zoom Out');             
        }
    };

    $cleanUpMenu->($menu_item, $self);

    foreach my $name (list_xclient_names()){
        $addSubMenuItem->($self, $menu_item, $name);
    }

    return ;
}


#===========================================================

sub RECEIVE_FILTER{
    my ($_connect, $_request, $_obj, @list) = @_;
    my ($_status, $_response) 
        = (404, $_obj->zMapZmapConnector->basic_error("Unknown Command"));

    my $lookup = {register_client => 'zMapRegisterClient'};
    my $action = '';
    my $reqXML = parse_request($_request);
    $action    = $reqXML->{'action'};

    foreach my $valid(@list){
        if($action eq $valid
           && ($valid = $lookup->{$valid}) # N.B. THIS SHOULD BE ASSIGNMENT NOT EQUALITY 
           && $_obj->can($valid)){
            ($_status, $_response) 
                = $_obj->$valid($_request);
            last;
        }
    }
    return ($_status, $_response) ;
}

sub open_clones{
    my ($watch) = @_;
    my ($self)  = @{$watch->Args('-store')};

    my @clones = $self->clone_list;
    my ($chr, $st, $end) = split(/\.|\-/, $clones[0]);

    my $seg = newXMLObj('segment');
    setObjNameValue($seg, 'sequence', "@clones");
    setObjNameValue($seg, 'start', 1);
    setObjNameValue($seg, 'end', '0');

    $self->zMapMakeRequest($seg, 'new');

    $self->zMapUpdateMenu();
    
    $watch->Unwatch;
}

sub RESPONSE_HANDLER{
    my ($self, $xml) = @_;
    warn "In RESPONSE_HANDLER\n" if $ZMAP_DEBUG;
    my ($name, $id) = ($xml->{'response'}->{'zmapid'}, $xml->{'response'}->{'windowid'});
    if($name){
        xclient_with_name($name, $id) if $id;
        $self->zMapSetEntryValue($name);
    }
    return ;
}
sub ERROR_HANDLER{
    my ($self, $status, $xml) = @_;
    $xml = $xml->{'error'}; # this is all we care about
    warn "In ERROR_HANDLER\n" if $ZMAP_DEBUG;
    if($status == 400){

    }elsif($status == 401){

    }elsif($status == 402){

    }elsif($status == 403){

    }elsif($status == 404){
        # could do something clever here so that we don't send the same window this command again.
        my $message = $xml->{'message'};
    }elsif($status == 412){
        delete_xclient_with_id($xml->{'windowid'});
    }elsif($status == 500){
    
    }elsif($status == 501){
    
    }elsif($status == 502){
    
    }elsif($status == 503){
    
    }else{
        warn "I know nothing about status $status\n";
    }
    return ;
}



1;

__END__


=pod

=head1 REMOVE

sub MenuCanvasWindow::XaceSeqChooser::update_display{
    my ($self , $ace) = @_ ;
    warn "cannot update display yet";
    return 1;
}
# this should be called when a user tries to save, but no ZMap is opened
sub MenuCanvasWindow::XaceSeqChooser::open_dialogue{
    my ($self) = @_ ;
    warn "cannot open dialogue yet";
    return undef;

}

sub MenuCanvasWindow::XaceSeqChooser::get_window_id {
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


sub MenuCanvasWindow::XaceSeqChooser::attach {
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


=head1 NAME - MenuCanvasWindow::ZmapSeqChooser

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk


=cut

__DATA__



