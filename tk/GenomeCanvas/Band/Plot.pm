
### GenomeCanvas::Band::Plot

package GenomeCanvas::Band::Plot;

use strict;
use Carp;
use GenomeCanvas::Band;
use Symbol 'gensym';

use vars '@ISA';
@ISA = ('GenomeCanvas::Band');

sub plot_method {
    my( $band, $plot_method ) = @_;
    
    if ($plot_method) {
        $band->{'_plot_method'} = $plot_method;
    }
    return $band->{'_plot_method'} || "gc_profile";
}

sub height {
    my( $band, $height ) = @_;
    
    if ($height) {
        $band->{'_height'} = $height;
    }
    return $band->{'_height'} || $band->font_size * 10;
}

sub range {
    my( $band, $low, $high ) = @_;
    
    if (defined $low) {
        confess "Invalid arguments"
            unless defined($high) and $low < $high;
        $band->{'_range'} = [$low,$high];
    }
    if (my $range = $band->{'_range'}) {
        return @$range;
    } else {
        return (0,1);
    }
}


sub show_midline {
    my( $self, $show_midline ) = @_;
    
    if (defined $show_midline) {
        $self->{'_show_midline'} = $show_midline ? 1 : 0;
    }
    $show_midline = $self->{'_show_midline'};
    if (defined $show_midline) {
        return $show_midline;
    } else {
        return 1;
    }
}


sub render {
    my( $band ) = @_;
    
    # Draw grey boxes in plot where there is sequence missing
    $band->draw_sequence_gaps;
    
    # Draw the plot where there is sequence
    $band->create_plot;
    
    #$band->draw_cpg_islands;
    
    # Draw the axes on top of the plot
    $band->draw_plot_axes;
}


sub draw_cpg_islands {
    my( $band ) = @_;
    
    my $vc = $band->virtual_contig
        or confess "No virtual contig attached";

    my $height    = $band->height;
    my $canvas    = $band->canvas;
    my $y_offset  = $band->y_offset;
    my $rpp       = $band->residues_per_pixel;
    my @tags      = $band->tags;

    my $half_height = $height / 2;
    my $y_middle = $y_offset + $half_height;
    
    my $color = '#4169e7';
    #my $max_score = log(5000);
    my $max_score = 5000;
    foreach my $cpg ($vc->get_all_SimpleFeatures_by_feature_type('cpg_island')) {
        my $x1 = $cpg->start / $rpp;
        my $x2 = $cpg->end   / $rpp;
        #my $score = log($cpg->score);
        my $score = $cpg->score;
        my $cpg_height = $half_height * ($score / $max_score);
        
        my $y1 = $y_middle - $cpg_height;
        my $y2 = $y_middle + $cpg_height;
        
        $canvas->createRectangle(
            $x1, $y1, $x2, $y2,
            -fill       => $color,
            -outline    => $color,
            -width      => 0,
            -tags       => [@tags],
            );
    }
}


sub band_color {
    my( $self, $band_color ) = @_;
    
    if ($band_color) {
        $self->{'_band_color'} = $band_color;
    }
    return $self->{'_band_color'} || 'red';
}


sub create_plot {
    my( $band ) = @_;
    
    my $plot_method = $band->plot_method;
    
    my @tags = $band->tags;
    my $rpp = $band->residues_per_pixel;
    my $vc = $band->virtual_contig
        or confess "No virtual contig attached";
    my $seq_length = $vc->length;
    my $canvas = $band->canvas;
    my $height = $band->height;
    my $y_offset = $band->y_offset;
    my $y_max = $y_offset + $height;
    my $y_middle = $y_offset + ($height / 2);
    
    
    # Get the map of the chunks of sequence
    my @seq_coord = $band->sequence_chunk_coords;

    my( $low, $high ) = $band->range;
    my $scale = $height / ($high - $low);
    my $fill = $band->band_color;
    for (my $i = 0; $i < @seq_coord; $i += 2) {
        my( $start, $end ) = @seq_coord[$i, $i+1];

        my( @plot );
        my( $pos, $value ) = $band->$plot_method($start, $end);
        for (my $j = 0; $j < @$pos; $j++) {
            my $p = $pos->[$j];
            my $x = $p / $rpp;

            my $v = $value->[$j];
            my $y = $y_max - (($v - $low) * $scale);
            push(@plot, $x, $y);
        }

        # Draw midline
        if ($band->show_midline) {
            my $x1 = $start / $rpp;
            my $x2 = $end   / $rpp;
            $canvas->createLine(
                $x1, $y_middle, $x2, $y_middle,
                #-fill       => '#ff6666',
                -fill       => '#666666',
                -width      => 1,
                -tags       => [@tags],
                );
        }

        # Draw plot
        $canvas->createLine(
            @plot,
            -fill       => $fill,
            -width      => 2,
            -tags       => [@tags],
            );
    }
}


sub pos_value_file {
    my( $self, $pos_value_file ) = @_;
    
    if ($pos_value_file) {
        $self->{'_pos_value_file'} = $pos_value_file;
    }
    return $self->{'_pos_value_file'};
}

sub pos_value {
    my( $self, $start, $end ) = @_;
    
    my $file = $self->pos_value_file
        or confess "pos_value_file not set";
    my $fh = gensym();
    open $fh, $file or confess "Can't read '$file' : $!";
    
    my $pos_list   = [];
    my $value_list = [];
    if ($start and $end) {
        while (<$fh>) {
            my ($pos, $value) = split;
            next if $pos < $start or $pos > $end;
            push(@$pos_list, $pos);
            push(@$value_list, $value);
        }
    } else {
        while (<$fh>) {
            my ($pos, $value) = split;
            push(@$pos_list, $pos);
            push(@$value_list, $value);
        }
    }
    return($pos_list, $value_list);
}

sub gc_profile {
    my( $band, $start, $end ) = @_;
    
    my $rpp = $band->residues_per_pixel;
    
    my $vc = $band->virtual_contig;
    # Get the sequence for this chunk, and
    # check that it's the correct length.
    my $seq_length = $end - $start + 1;
    my $seq_str = $vc->subseq($start, $end);
    unless ($seq_length == length($seq_str)) {
        confess "subseq string is ", length($seq_str), " long, not $seq_length";
    }
    $seq_str = lc $seq_str;
    
    my( @pos, @value );
    my $tile_incr = $rpp * 2;
    my $tile_length = $tile_incr * 2;
    for (my $j = 0; $j < $seq_length; $j += $tile_incr) {
        my $tile_seq = substr($seq_str, $j, $tile_length);
        my $length = length($tile_seq);
        if ($length < $tile_length) {
            # Probably at the end of the sequence
            my $last_tile_seq = substr($seq_str, -1 * $tile_length, $tile_length);
            $length = length($last_tile_seq);
            if ($length == $tile_length) {
                my $gc_count = $last_tile_seq =~ tr/gc/gc/;
                my $gc_fraction = $gc_count / $length;
                push(@pos, $start + $seq_length - $tile_incr);
                push(@value, $gc_fraction);
            } else {
                warn "Didn't get full tile ($length != $tile_length)";
            }
        } else {
            my $gc_count = $tile_seq =~ tr/gc/gc/;
            my $gc_fraction = $gc_count / $length;
            push(@pos, $start + $j + $tile_incr);
            push(@value, $gc_fraction);
        }
    }
    
    return (\@pos, \@value);
}

sub x_major {
    my( $band, $x_major ) = @_;
    
    if ($x_major) {
        $band->{'_x_major'} = $x_major;
    }
    return $band->{'_x_major'} || 0.5;
}

sub draw_plot_axes {
    my( $band ) = @_;
    
    my $y_offset    = $band->y_offset;
    my $height      = $band->height;
    my $width       = $band->width;
    my @tags        = $band->tags;
    
    my $y_max = $y_offset + $height;
    my( $low, $high ) = $band->range;
    my $scale = $height / ($high - $low);
    my $major = $band->x_major;

    # Y axis ticks
    for (my $i = $low; $i <= $high; $i += $major) {
        $i = sprintf("%.1f", $i);
        my $y = $y_max - (($i - $low) * $scale);
        # Left axis
        $band->tick_label($i, 'w', 0,      $y);
        # Right axis
        $band->tick_label($i, 'e', $width, $y);
    }
    
    # Scale along bottom of plot
    my $vc = $band->virtual_contig;
    my $chr_start  = $vc->_global_start;
    my $seq_length = $vc->length;
    my $seq_end = $chr_start + $seq_length;
    my $rpp = $band->residues_per_pixel;
    # Choose minumum space between labeled ticks based on font size
    my $min_interval = $rpp * $band->font_size * 4;
    my $interval = 1;
    my $multiplier = 5;
    while ($interval < $min_interval) {
        $interval *= $multiplier;
        $multiplier = ($multiplier == 5) ? 2 : 5;
    }
    my $precision = 7 - length($interval);
    #warn "interval = $interval, precision = $precision";
    {
        my( $i );
        if (my $rem = $chr_start % $interval) {
            $i = $chr_start - $rem + $interval;
        } else {
            $i = $chr_start;
        }
        #warn "First label = $i";
        for (; $i <= $seq_end; $i += $interval) {
            my $Mbp = sprintf("%.${precision}f", $i / 1e6);
            $band->tick_label($Mbp, 's', ($i - $chr_start) / $rpp, $y_max);
        }
    }

    # Rectangle in front of everything else
    my $canvas = $band->canvas;
    my $outline = $canvas->createRectangle(
        0, $y_offset, $width, $y_max,
        -fill       => undef,
        #-border     => 'black',
        -width      => 1,
        -tags       => [@tags],
        );
}

sub sequence_chunk_coords {
    my( $band ) = @_;
    
    my $vc = $band->virtual_contig;
    # Make a map of the chunks of sequence
    my( @seq_coord );
    foreach my $map_c ($vc->_vmap->each_MapContig) {
        my $start  = $map_c->start;
        my $end    = $map_c->end;

        if (@seq_coord) {
            if ($seq_coord[$#seq_coord] == $start - 1) {
                $seq_coord[$#seq_coord] = $end;
            } else {
                push(@seq_coord, $start, $end);
            }
        } else {
            @seq_coord = ($start, $end);
        }
    }
    return @seq_coord;
}

1;

__END__

=head1 NAME - GenomeCanvas::Band::Plot

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

