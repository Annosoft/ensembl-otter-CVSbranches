
### GenomeCanvas::Band::RegionFeatures

package GenomeCanvas::Band::RegionFeatures;

use strict;
use strict;
use Carp;
use base 'GenomeCanvas::Band';


sub new {
    my( $pkg, @args ) = @_;
    
    my $self = $pkg->SUPER::new(@args);
    $self->height(60);
    return $self;
}

sub analysis_type {
    my( $self, $analysis_type ) = @_;
    
    if ($analysis_type) {
        $self->{'_analysis_type'} = $analysis_type;
    }
    return $self->{'_analysis_type'} || 'stats_region';
}

sub get_label_start_end {
    my( $self ) = @_;
    
    my $vc = $self->virtual_contig
        or confess "No virtual contig attached";
    
    my $type = $self->analysis_type;
    my $feat_startend = {};
    foreach my $feat ($vc->get_all_SimpleFeatures_by_feature_type($type)) {
        my $name = $feat->source_tag;
        #print STDERR "$name\n";
        if ($name =~ s/_left$//) {
            my $a = $feat_startend->{$name} ||= [];
            $a->[0] = $feat->start;
        }
        elsif ($name =~ s/_right$//) {
            my $a = $feat_startend->{$name} ||= [];
            $a->[1] = $feat->end;
        }
    }
    return $feat_startend;
}

sub render {
    my( $self ) = @_;
    
    my $vc = $self->virtual_contig
        or confess "No virtual contig attached";
    my $global_offset = $vc->_global_start - 1;
    confess "Can't get global_start from Virtual Contig"
        unless defined($global_offset);

    my $height    = $self->height;
    my $canvas    = $self->canvas;
    my $y_offset  = $self->y_offset;
    my $rpp       = $self->residues_per_pixel;
    my $color     = $self->band_color;
    my $pad       = $self->font_size / 2;
    my @tags      = $self->tags;
    
    $canvas->createRectangle(
        0, $y_offset, $self->width, $y_offset + $height,
        -fill       => undef,
        -outline    => undef,
        -tags       => [@tags],
        );

    my $lse = $self->get_label_start_end;

    my $y1 = $y_offset;
    my $y2 = $y_offset + $height;
    my $text_y = $y_offset + (0.6 * $height);
    
    my @even_odd_colors = (
        {
            'fill' => '#999999',
            'text' => 'black',
        },
        {
            'fill' => 'black',
            'text' => 'white',
        },
        );
    my @region_names = sort keys %$lse;
    my $font = ['helvetica', 0.8 * $height, 'bold'];
    for (my $i = 0; $i < @region_names; $i++) {
        my $region = $region_names[$i];
        my ($text) = $region =~ /(\d+)/;
        my ($start, $end) = @{$lse->{$region}};
        my $color = $even_odd_colors[$i % 2];
        
        my $x1 = ($start - $global_offset) / $rpp;
        my $x2 = ($end   - $global_offset) / $rpp;
        my $text_x = $x1 + (($x2 - $x1) / 2);
        
        $canvas->createRectangle(
            $x1, $y_offset, $x2, $y2,
            -fill       => $color->{'fill'},
            -outline    => undef,
            -tags       => [@tags],
            );

        foreach my $x_pos ($x1, $x2) {
            $canvas->createLine(
                $x_pos, $y_offset - $pad, $x_pos, $y2 + $pad,
                -fill       => 'black',
                -width      => 2,
                -tags       => [@tags],
                );
        }
        
        $canvas->createText(
            $text_x, $text_y,
            -anchor     => 'center',
            -font       => $font,
            -text       => $text,
            -fill       => $color->{'text'},
            -tags       => [@tags],
            );
    }
    
    #while (<$fh>) {
    #    # Cunning or what!
    #    my ($start, $end) = map $_ -= $global_offset, (split)[@ind];
    #    warn "Start End = $start\t$end\n";
    #
    #    my $x1 = $start / $rpp;
    #    my $x2 = $end   / $rpp;
    #    
    #    $canvas->createRectangle(
    #        $x1, $y_offset, $x2, $y2,
    #        -fill       => $color,
    #        -outline    => $color,
    #        -width      => 0.5,
    #        -tags       => [@tags],
    #        );
    #}
    
}



1;

__END__

=head1 NAME - GenomeCanvas::Band::RegionFeatures

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

