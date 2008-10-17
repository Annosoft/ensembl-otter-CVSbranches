
### EditWindow::Exonerate

package EditWindow::Exonerate;

use strict;
use warnings;
use Carp;

use Hum::Pfetch;
use Hum::FastaFileIO;
use Bio::Otter::Lace::Exonerate;
use Tk::LabFrame;
use Tk::FileDialog;
use Tk::Balloon;

use base 'EditWindow';

my $PROT_SCORE = 150;
my $DNA_SCORE  = 2000;
my $DNAHSP     = 120;

sub initialise {
	my ( $self ) = @_;
	my $top  = $self->top;
	my $exon = Bio::Otter::Lace::Exonerate->new;
	$exon->AceDatabase($self->XaceSeqChooser->AceDatabase);
	$self->exonerate($exon);
	### Query frame
	my $query_frame = $top->LabFrame(
									  -borderwidth => 3,
									  -label       => 'Query sequences',
									  -labelside   => 'acrosstop',
	)->pack( -side => 'top', );
	## Molecule type
	my $type_frame =
	  $query_frame->Frame( -border => 3, )->pack( -side => 'top', );
	$type_frame->Label(
						-text   => 'Molecule type:',
						-anchor => 's',
						-padx   => 6,
	)->pack( -side => 'left' );

	# Pad
	$type_frame->Frame( -width => 20, )->pack( -side => 'left' );

	# dna or protein sequence
	my $type = 'dna';
	$exon->query_type($type);
	my $mol_type = sub {
		if ( $type eq 'dna' ) {
			$self->set_entry( 'score', $DNA_SCORE );
			$self->dnahsp->configure( -state => 'normal' );
			$self->set_entry( 'dnahsp', $DNAHSP );
			$self->set_entry( 'method_tag', 'Exon_DNA' );
			$self->set_entry( 'logic_name', 'Exon_DNA' );
		}
		else {
			$self->set_entry( 'score',  $PROT_SCORE );
			$self->set_entry( 'dnahsp', 0 );
			$self->dnahsp->configure( -state => 'disable' );
			$self->set_entry( 'method_tag', 'Exon_PROT' );
			$self->set_entry( 'logic_name', 'Exon_PROT' );
		}
		$exon->query_type($type);
	};
	foreach (qw/dna protein/) {
		$type_frame->Radiobutton(
								  -text     => $_,
								  -variable => \$type,
								  -value    => $_,
								  -command  => $mol_type
		)->pack( -side => 'right' );
	}
	## Accession entry box
	my $match_frame =
	  $query_frame->Frame( -border => 3, )->pack( -side => 'top', );
	$match_frame->Label(
						 -text   => 'Accessions:',
						 -anchor => 's',
						 -padx   => 6,
	)->pack( -side => 'left' );
	$self->match(
				$match_frame->Entry( -width => 24, )->pack( -side => 'left' ) );
	## Fasta file entry box
	my $fname;
	my $Horiz      = 1;
	my $file_frame =
	  $query_frame->Frame( -border => 3, )->pack( -side => 'top', );
	$file_frame->Label(
						-text   => 'OR Fasta file:',
						-anchor => 's',
						-padx   => 6,
	)->pack( -side => 'left' );

	# Pad between entries
	$file_frame->Frame( -width => 10, )->pack( -side => 'top' );
	my $LoadDialog = $file_frame->FileDialog(
											  -Title   => 'Select a Fasta file',
											  -Create  => 0,
											  -FPat    => '*fa',
											  -ShowAll => 'NO'
	);
	$self->fasta_file(
		 $file_frame->Entry( -textvariable => \$fname )->pack( -side => 'left' ) );
	$file_frame->Button(
		-text    => 'Browse...',
		-command => sub {
			$fname = $LoadDialog->Show( -Horiz => $Horiz );
		}
	)->pack( -side => 'left' );
	## Sequence text box
	my $txt_frame =
	  $query_frame->Frame( -border => 3, )->pack( -side => 'top', );
	$txt_frame->Label(
					   -text   => 'OR Fasta sequences',
					   -anchor => 's',
					   -padx   => 6,
	)->pack( -side => 'top' );
	$self->fasta_txt(
		   $txt_frame->Scrolled( "Text", -background => 'white', -height => 12 )
			 ->pack( -side => 'top', ) );
	### Exonerate and diplay parameters
	my $param_frame = $top->LabFrame(
									  -borderwidth => 3,
									  -label       => 'Parameters',
									  -labelside   => 'acrosstop',
	)->pack( -side => 'top', );
	## Score and dna hsp thresholds
	my $threshold_frame =
	  $param_frame->Frame( -border => 3, )->pack( -side => 'top', );
	$threshold_frame->Label(
							 -text   => 'Threshold ',
							 -anchor => 's',
							 -padx   => 6,
	)->pack( -side => 'left' );
	$threshold_frame->Label(
							 -text   => 'Score:',
							 -anchor => 's',
							 -padx   => 6,
	)->pack( -side => 'left' );
	$self->score(
				  $threshold_frame->Entry(
										   -width   => 9,
										   -justify => 'right',
					)->pack( -side => 'left' )
	);
	$self->set_entry( 'score', $DNA_SCORE );
	$threshold_frame->Label(
							 -text   => 'Dna hsp:',
							 -anchor => 's',
							 -padx   => 6,
	)->pack( -side => 'left' );
	$self->dnahsp(
				   $threshold_frame->Entry(
											-width   => 9,
											-justify => 'right',
					 )->pack( -side => 'left' )
	);
	$self->set_entry( 'dnahsp', $DNAHSP );
	## ace method tag/color and logic_name
	my $display_frame =
	  $param_frame->Frame( -border => 3, )->pack( -side => 'top', );
	$display_frame->Label(
						   -text   => 'Method ',
						   -anchor => 's',
						   -padx   => 6,
	)->pack( -side => 'left' );
	$display_frame->Label(
						   -text   => 'Tag:',
						   -anchor => 's',
						   -padx   => 6,
	)->pack( -side => 'left' );
	$self->method_tag(
					   $display_frame->Entry(
											  -width   => 9,
											  -justify => 'right',
						 )->pack( -side => 'left' )
	);
	$display_frame->Label(
						   -text   => 'Color:',
						   -anchor => 's',
						   -padx   => 6,
	)->pack( -side => 'left' );

	$self->method_color('red');
	my $button_color;
	my $balloon    = $display_frame->Balloon(-bg => "khaki2", -initwait => 500);
	my $show = sub { $self->show_color_panel($button_color); };
	$button_color = $display_frame->Button(
						   -background => 'red',
						   -activebackground => 'red',
						   -relief => 'flat',
						   -borderwidth => 1,
						   -width   => 5,
						   -command => $show
	)->pack( -side => 'left' );
	$balloon->attach($button_color, -msg => "click to change color");


	$display_frame->Label(
						   -text   => 'Logic_name:',
						   -anchor => 's',
						   -padx   => 6,
	)->pack( -side => 'left' );
	$self->logic_name(
					   $display_frame->Entry(
											  -width   => 9,
											  -justify => 'right',
						 )->pack( -side => 'left' )
	);
	$self->set_entry( 'method_tag',   'Exon_DNA' );
	$self->set_entry( 'logic_name',   'Exon_DNA' );
	### Commands
	my $button_frame = $top->Frame->pack(    -side => 'top',
										  -fill => 'x', );
	my $launch = sub {
		$self->launch_exonerate or return;
		$top->withdraw;
	};
	$button_frame->Button(
						   -text      => 'Launch',
						   -underline => 0,
						   -command   => $launch,
	)->pack( -side => 'left' );
	$top->bind( '<Control-l>', $launch );
	$top->bind( '<Control-L>', $launch );

	# Update coords
	my $update = sub {
		$self->update_from_clipboard;
	};
	$button_frame->Button(
						   -text      => 'Update',
						   -underline => 0,
						   -command   => $update,
	)->pack( -side => 'left' );
	$top->bind( '<Control-u>', $update );
	$top->bind( '<Control-U>', $update );

	# Manage window closes and destroys
	my $close_window = sub { $top->withdraw };
	$button_frame->Button(
						   -text    => 'Close',
						   -command => $close_window,
	)->pack( -side => 'right' );
	$top->bind( '<Control-w>', $close_window );
	$top->bind( '<Control-W>', $close_window );
	$top->protocol( 'WM_DELETE_WINDOW', $close_window );
	$top->bind( '<Destroy>', sub {	$self = undef ;  } );
}

sub update_from_XaceSeqChooser {
	my ( $self ) = @_;
	$self->update_from_clipboard;
	my $top = $self->top;
	$top->deiconify;
	$top->raise;
}

sub query_Sequence {
	my ( $self, $query_Sequence ) = @_;
	if ($query_Sequence) {
		$self->{'_query_Sequence'} = $query_Sequence;
	}
	return $self->{'_query_Sequence'};
}

sub update_from_clipboard {
	my ($self) = @_;
	if ( my ( $name, $start, $end ) = $self->name_start_end_from_fMap_blue_box )
	{
		$self->set_entry( 'match', $name );
	}
}

sub set_entry {
	my ( $self, $method, $txt ) = @_;
	my $entry = $self->$method();
	$entry->delete( 0, 'end' );
	$entry->insert( 0, $txt );
}

sub get_entry {
	my ( $self, $method ) = @_;
	my $txt = $self->$method()->get or return;
	$txt =~ s/\s//g;
	return $txt;
}

sub method_tag {
	my ( $self, $match ) = @_;
	if ($match) {
		$self->{'_method_tag'} = $match;
	}
	return $self->{'_method_tag'};
}

sub fasta_txt {
	my ( $self, $txt ) = @_;
	if ($txt) {
		$self->{'_fasta_txt'} = $txt;
	}
	return $self->{'_fasta_txt'};
}

sub fasta_file {
	my ( $self, $file ) = @_;
	if ($file) {
		$self->{'_fasta_file'} = $file;
	}
	return $self->{'_fasta_file'};
}

sub method_color {
	my ( $self, $match ) = @_;
	if ($match) {
		$self->{'_method_color'} = $match;
	}
	return $self->{'_method_color'};
}

sub logic_name {
	my ( $self, $match ) = @_;
	if ($match) {
		$self->{'_logic_name'} = $match;
	}
	return $self->{'_logic_name'};
}

sub score {
	my ( $self, $match ) = @_;
	if ($match) {
		$self->{'_score'} = $match;
	}
	return $self->{'_score'};
}

sub dnahsp {
	my ( $self, $match ) = @_;
	if ($match) {
		$self->{'_dnahsp'} = $match;
	}
	return $self->{'_dnahsp'};
}

sub match {
	my ( $self, $match ) = @_;
	if ($match) {
		$self->{'_match'} = $match;
	}
	return $self->{'_match'};
}

sub genomic {
	my ( $self, $genomic ) = @_;
	if ($genomic) {
		$self->{'_genomic'} = $genomic;
	}
	return $self->{'_genomic'};
}

sub genomic_start {
	my ( $self, $genomic_start ) = @_;
	if ( defined $genomic_start ) {
		$self->{'_genomic_start'} = $genomic_start;
	}
	return $self->{'_genomic_start'};
}

sub genomic_end {
	my ( $self, $genomic_end ) = @_;
	if ( defined $genomic_end ) {
		$self->{'_genomic_end'} = $genomic_end;
	}
	return $self->{'_genomic_end'};
}

sub flank {
	my ( $self, $flank ) = @_;
	if ($flank) {
		$self->{'_flank'} = $flank;
	}
	return $self->{'_flank'};
}

sub revcomp_ref {
	my ( $self, $revcomp_ref ) = @_;
	if ($revcomp_ref) {
		$self->{'_revcomp_ref'} = $revcomp_ref;
	}
	return $self->{'_revcomp_ref'};
}

sub exonerate {
	my ( $self, $exonerate ) = @_;
	if ($exonerate) {
		$self->{'_exonerate'} = $exonerate;
	}
	return $self->{'_exonerate'};
}

sub show_color_panel {
	my ($self,$button_color) = @_;
	my $color_dialog = $self->top()->DialogBox( -title => 'Choose a color', -buttons => ['Cancel'] );
	my @color_name = Hum::Ace::Colors::list_all_color_names_by_value();
	my $balloon    = $color_dialog->Balloon(-bg => "khaki2", -initwait => 700);
	for (my $i = 0; $i < @color_name; $i++) {
       my $name = $color_name[$i];
        my $hex = Hum::Ace::Colors::acename_to_webhex($name);
        my $color = sub {
        	$button_color->configure(
        					-background => $hex,
        					-activebackground => $hex);
			$self->method_color($name);
			$color_dialog->Exit();
        };
        my $button = $color_dialog->add("Button",
        				   -background => $hex,
						   -activebackground => $hex,
						   -relief => 'flat',
						   -borderwidth => 1,
						   -width   => 3,
						   -command => $color)->pack;
		$balloon->attach($button, -msg => "$name");


	}
	$color_dialog->bind( '<Destroy>', sub { $self = undef } );

	return $color_dialog->Show();
}

sub XaceSeqChooser {
	my ($self,$xc) = @_;
	if ($xc) {
		$self->{'_xc'} = $xc;
	}
	return $self->{'_xc'};
}

sub launch_exonerate {
	my ($self) = @_;
	my $seq = $self->get_query_seq();
	print STDOUT "Found " . scalar(@$seq) . " sequences\n";
	my $score   = $self->get_entry('score');
	my $dnahsp  = $self->get_entry('dnahsp');
	my $m_tag   = $self->get_entry('method_tag');
	my $m_color = $self->method_color();
	my $l_name  = $self->get_entry('logic_name');
	unless ( $score and $m_tag and $m_color and $l_name and $seq) {
		warn "Missing parameters\n";
		return;
	}
	my $exonerate = $self->exonerate;
	$exonerate->query_seq($seq);
	$exonerate->score($score);
	$exonerate->dnahsp($dnahsp);
	$exonerate->method_tag($m_tag);
	$exonerate->method_color($m_color);
	$exonerate->logic_name($l_name);
	my $seq_file = $exonerate->write_seq_file();
	if($seq_file){
		$exonerate->initialise($seq_file);
		my $ace_text = $exonerate->run or return;
		# delete query file
		unlink $seq_file;
		$self->top->Busy;
		# Need to add new method to collection if we don't have it already
    	my $coll = $exonerate->AceDatabase->MethodCollection;
    	my $coll_zmap = $self->XaceSeqChooser->Assembly->MethodCollection;
    	my $method = $exonerate->ace_Method;
    	unless ($coll->get_Method_by_name($method->name) ||
    			$coll_zmap->get_Method_by_name($method->name)) {
        	$coll->add_Method($method);
        	$coll_zmap->add_Method($method);
        	$self->XaceSeqChooser->save_ace($coll->ace_string());
    	}

		$self->XaceSeqChooser->save_ace($ace_text);
		$self->XaceSeqChooser->zMapWriteDotZmap;
		$self->XaceSeqChooser->resync_with_db();
		$self->XaceSeqChooser->zMapLaunchZmap;
		$self->top->Unbusy;
	}
}

my $seq_tag = 1;

sub get_query_seq {
	my ($self) = @_;
	my @seq;
	if ( $self->get_entry('match') ) {
		my @accessions = split /\,|\;/, $self->get_entry('match');
		if (@accessions) {
			push @seq, Hum::Pfetch::get_Sequences(@accessions);
		}
	}
	if ( my $string = $self->fasta_txt->get( '1.0', 'end' ) ) {
		if( $string =~ /\S/ && !($string =~ />/)) {
			$string = ">Unknow_$seq_tag\n".$string; $seq_tag++;
			push @seq, Hum::FastaFileIO->new_String_IO($string)->read_all_sequences;
		}
	}
	if ( $self->get_entry('fasta_file') ) {
		push @seq, Hum::FastaFileIO->new( $self->get_entry('fasta_file') )
		  ->read_all_sequences;
	}
	return \@seq;
}

sub name_start_end_from_fMap_blue_box {
	my ($self) = @_;
	my $tk = $self->top;
	my $text = $self->get_clipboard_text or return;

	#warn "clipboard: $text";
	# Match fMap "blue box"
	if ( $text =~
/^(?:<?(?:Protein|Sequence)[:>]?)?\"?([^\"\s]+)\"?\s+-?(\d+)\s+-?(\d+)\s+\(\d+\)/
	  )
	{
		my $name  = $1;
		my $start = $2;
		my $end   = $3;
		( $start, $end ) = ( $end, $start ) if $start > $end;

		#warn "Got ($name, $start, $end)";
		return ( $name, $start, $end );
	}
	else {
		return;
	}
}

sub DESTROY {
    my ($self) = shift;

	warn "Freeing exonerateWindow '$self'\n";
}

1;

__END__

=head1 NAME - EditWindow::Exonerate

=head1 AUTHOR

Anacode B<email> anacode@sanger.ac.uk
