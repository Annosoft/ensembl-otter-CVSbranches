
### EditWindow::LocusName

package EditWindow::LocusName;

use strict;
use warnings;
use Carp;

use Tk::SmartOptionmenu;
use base 'EditWindow';
use Scalar::Util 'weaken';


sub initialise {
    my ($self) = @_;
    
    my $top = $self->top;

    my $name_frame = $top->Frame(
        -border     => 3,
        )->pack(-side => 'top', -fill => 'x');

    # Menu for chosing from list of existing locus names
    my ($chosen, $menu_list) = $self->make_menu_choices;
    $name_frame->Label(-text => 'Rename locus:')->pack(-side => 'left', -padx => 3);
    my $name_menu = $name_frame->SmartOptionmenu(
        -options    => $menu_list,
        -variable   => \$chosen,
    )->pack(-side => 'left');
    $self->{'_chosen_name'} = \$chosen;
    
    $name_frame->Label(-text => 'to:')->pack(-side => 'left', -padx => 3);
    
    # Entry for editing new name of loucs
    my $new_name = $name_frame->Entry(-width => 20)->pack(-side => 'left');
    $self->{'_new_name_entry'} = $new_name;

    # Command which copies the chosen name to the new name Entry widget
    my $put_chosen_in_entry = sub {
        $new_name->delete(0, 'end');
        $new_name->insert(0, ${$self->{'_chosen_name'}});
        $new_name->selectionRange(0, 'end');
        $new_name->focus;
    };
    $name_menu->configure(-command => $put_chosen_in_entry);
    $put_chosen_in_entry->();
    
    my $button_frame = $top->Frame(
        -border     => 3,
        )->pack(-side => 'top', -fill => 'x');
    
    # Button which renames the locus
    my $do_rename = sub { $self->do_rename };
    my $rename_button = $button_frame->Button(
            -text       => 'Rename',
            -default    => 'active',
            -command    => $do_rename,
            )->pack(-side => 'left');
    my $press_rename_button = sub{
        $rename_button->focus;
        $rename_button->invoke;
        };
    $top->bind('<Return>',      $press_rename_button);
    $top->bind('<KP_Enter>',    $press_rename_button);
    
    my $cancel = sub { $top->destroy };
    $button_frame->Button(-text => 'Cancel', -command => $cancel)->pack(-side => 'right');
    $top->bind('<Escape>', $cancel);
    

    
    $top->bind('<Destroy>', sub{ $self = undef; });

    return;
}

sub do_rename {
    my ($self) = @_;
    
    my $old_name = $self->chosen_name;
    my $new_name = $self->get_new_name;
    if ($old_name eq $new_name) {
        return;
    }
    warn "Renaming Locus '$old_name' to '$new_name'\n";
    
    $self->XaceSeqChooser->do_rename_locus($old_name, $new_name);

    $self->top->destroy;

    return;
}

sub locus_name_arg {
    my ($self, $arg) = @_;
    
    if ($arg) {
        $self->{'_locus_name_arg'} = $arg;
    }
    return $self->{'_locus_name_arg'};
}

sub chosen_name {
    my ($self) = @_;
    
    return ${$self->{'_chosen_name'}};
}

sub get_new_name {
    my ($self) = @_;
    
    return $self->{'_new_name_entry'}->get;
}

sub XaceSeqChooser {
    my ($self, $xc) = @_;
    
    if ($xc) {
        $self->{'_xace_seq_chooser'} = $xc;
        weaken($self->{'_xace_seq_chooser'});
    }
    return $self->{'_xace_seq_chooser'};
}

sub make_menu_choices {
    my ($self) = @_;
    
    my $xc = $self->XaceSeqChooser;
    my $sel_locus_name = $self->locus_name_arg;
    my @locus_name = $xc->list_Locus_names;

    # If we were passed a locus name, check that it is actually a locus name
    my $saw = 0;
    if ($sel_locus_name) {
        foreach my $name (@locus_name) {
            if ($name eq $sel_locus_name) {
                $saw = 1;
                last;
            }
        }
    }
    $sel_locus_name = undef unless $saw;
    
    # If we don't have a locus name, take the one from the first
    # selected subseq we find
    unless ($sel_locus_name) {
        foreach my $name ($xc->list_selected_subseq_names) {
            my $sub = $xc->get_SubSeq($name) or next;
            my $locus = $sub->Locus or next;
            $sel_locus_name = $locus->name;
            last;
        }
    }
    
    # Or choose the one at the top of the menu.
    $sel_locus_name ||= $locus_name[0];

    my $menu_list = [ map { [$_, $_] } @locus_name ];
    return ($sel_locus_name, $menu_list);
}

sub DESTROY {
    my ($self) = @_;
    
    warn "Destroying a '", ref($self), "'";

    return;
}

1;

__END__

=head1 NAME - EditWindow::LocusName

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

