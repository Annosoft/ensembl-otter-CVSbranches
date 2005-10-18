package BlastzAligner;

=head1 NAME

BlastzAligner.pm - module to do a whole genome alignment between two closely
related assemblies and create assembly entries from it.

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 AUTHOR

Patrick Meidl <pm2@sanger.ac.uk>

=head1 CONTACT

Post questions to the EnsEMBL development list ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use File::Path;

use constant FMT1 => "%-30s%10.0f (%3.2f%%)\n";
use constant FMT2 => "%-30s%10.0f\n";
use constant FMT3 => "%-8s%-12s%-5s%-10s%-10s%-10s%-10s\n";
use constant FMT4 => "%-8s%-12s%-5s%8.0f  %8.0f  %8.0f  %8.0f\n";

=head2 new

  Arg[-SUPPORT] : 
  Example       : 
  Description   : 
  Return type   : 
  Exceptions    : 
  Caller        : 

=cut

sub new {
    my ($class, @args) = @_;
    
    my $self = {};
    bless $self, $class;

    my ($support) = rearrange([qw(SUPPORT)], @args);
    $self->support($support);

    return $self;
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub create_tempdir {
    my $self = shift;

    # create tmpdir to store input and output
    my $user = `whoami`;
    chomp $user;
    my $tempdir = "/tmp/$user.".time.".".int(rand(1000));
    $self->support->log("Creating tmpdir $tempdir...\n");
    system("mkdir $tempdir") == 0 or
        $self->support->log_error("Can't create tmp dir $tempdir: $!\n");
    $self->support->log("Done.\n");

    $self->tempdir($tempdir);
    return $tempdir;
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub remove_tempdir {
    my $self = shift;
    rmtree($self->tempdir) or $self->support->log_warning("Could not delete ".$self->tempdir.": $!\n");
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub write_single_sequences {
    my ($self, $E_slice, $V_slice) = @_;

    unless (ref($E_slice) eq 'Bio::EnsEMBL::Slice' 
            and ref($V_slice) eq 'Bio::EnsEMBL::Slice') {
        $self->support->log_error("You must supply an Ensembl and a Vega Bio::EnsEMBL::Slice.\n");
    }
    
    my $tmpdir = $self->tempdir;
    my $bindir = $self->support->param('bindir');
    my $id = $self->id;

    my $E_fh = $self->support->filehandle('>', "$tmpdir/e_seq.$id.fa");
    my $V_fh = $self->support->filehandle('>', "$tmpdir/v_seq.$id.fa");

    print $E_fh join(':', ">e_seq.$id dna:chromfrag chromosome",
                          $self->support->param('ensemblassembly'),
                          $E_slice->start,
                          $E_slice->end,
                          $E_slice->strand
                    ), "\n";
    print $E_fh $E_slice->get_repeatmasked_seq(undef, 1)->seq, "\n";
    close($E_fh);

    print $V_fh join(':', ">v_seq.$id dna:chromfrag chromosome",
                          $self->support->param('assembly'),
                          $V_slice->start,
                          $V_slice->end,
                          $V_slice->strand
                    ), "\n";
    print $V_fh $V_slice->get_repeatmasked_seq(undef, 1)->seq, "\n";
    close($V_fh);

    # convert fasta to nib (needed for lavToAxt)
    system("$bindir/faToNib $tmpdir/e_seq.$id.fa $tmpdir/e_seq.$id.nib") == 0 or
        $self->support->log_error("Can't run faToNib: $!\n");
    system("$bindir/faToNib $tmpdir/v_seq.$id.fa $tmpdir/v_seq.$id.nib") == 0 or
        $self->support->log_error("Can't run faToNib: $!\n");
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub write_sequence {
    my ($self, $slice, $assembly, $basename1, $basename2) = @_;

    my $tmpdir = $self->tempdir;
    my $bindir = $self->support->param('bindir');

    my $fh = $self->support->filehandle('>', "$tmpdir/$basename1.fa");
    print $fh join(':', ">$basename1 dna:chromfrag chromosome",
                          $assembly,
                          $slice->start,
                          $slice->end,
                          $slice->strand
                    ), "\n";
    print $fh $slice->get_repeatmasked_seq(undef, 1)->seq, "\n";
    close($fh);
    
    # convert fasta to nib (needed for lavToAxt)
    system("$bindir/faToNib $tmpdir/$basename1.fa $tmpdir/$basename1.nib") == 0
        or $self->support->log_error("Can't run faToNib: $!\n");

    if ($basename2) {
        system("cat $tmpdir/$basename1.fa >> $tmpdir/$basename2.fa") == 0 or
            $self->support->log_error("Can't concat fasta files: $!\n");
    }
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub run_blastz {
    my ($self, $E_basename, $V_basename) = @_;

    my $tmpdir = $self->tempdir;
    my $bindir = $self->support->param('bindir');
    my $id = $self->id;

    my $blastz_cmd = qq($bindir/blastz $tmpdir/$E_basename.fa $tmpdir/$V_basename.fa Q=blastz_matrix.txt T=0 L=10000 H=2200 Y=3400 > $tmpdir/blastz.$id.lav);
    system($blastz_cmd) == 0 or
        $self->support->log_error("Can't run blastz: $!\n");
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub lav_to_axt {
    my $self = shift;

    my $tmpdir = $self->tempdir;
    my $bindir = $self->support->param('bindir');
    my $id = $self->id;

    system("$bindir/lavToAxt $tmpdir/blastz.$id.lav $tmpdir $tmpdir $tmpdir/blastz.$id.axt") == 0 or $self->support->log_error("Can't run lavToAxt: $!\n");
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub find_best_alignment {
    my $self = shift;

    my $tmpdir = $self->tempdir;
    my $bindir = $self->support->param('bindir');
    my $id = $self->id;

    system("$bindir/axtBest $tmpdir/blastz.$id.axt all $tmpdir/blastz.$id.best.axt") == 0 or $self->support->log_error("Can't run axtBest: $!\n");
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub parse_blastz_output {
    my $self = shift;

    # read file
    my $tmpdir = $self->tempdir;
    my $id = $self->id;
    my $fh = $self->support->filehandle('<', "$tmpdir/blastz.$id.best.axt");

    # initialize stats
    $self->init_stats(qw(match mismatch gap alignments bp));
    
    my $i = 1;
    my ($header, $e_seq, $v_seq);

    while (my $line = <$fh>) {
        # there are blocks of 4 lines, where line 1 is the header, line 2 is
        # e_seq, line3 is v_seq
        $header = $line unless (($i-1) % 4);
        $e_seq = $line unless (($i-2) % 4);
        chomp $e_seq;
        my @e_arr = split(//, $e_seq);
        $v_seq = $line unless (($i-3) % 4);
        chomp $v_seq;
        my @v_arr = split(//, $v_seq);

        # compare sequences letter by letter
        if ($i % 4 == 0) {
            my $match_flag = 0;
            $self->init_stats(qw(e_gap v_gap));
            my %coords;
            @coords{'v_id', 'e_start', 'e_end', 'v_start', 'v_end', 'strand'} =
                (split(/ /, $header))[4, 2, 3, 5, 6, 7];
            $coords{'v_id'} =~ s/v_seq\.(.*)/$1/;
            ($coords{'strand'} eq '+') ? ($coords{'strand'} = 1) :
                                         ($coords{'strand'} = -1);
            for (my $j = 0; $j < scalar(@e_arr); $j++) {
                # gap
                if ($e_arr[$j] eq '-' or $v_arr[$j] eq '-') {
                    $self->stats_incr('gap', 1);
                    $self->stats_incr('e_gap', 1) if ($e_arr[$j] eq '-');
                    $self->stats_incr('v_gap', 1) if ($v_arr[$j] eq '-');
                    $match_flag = 0;
                } else {
                    # match
                    if ($e_arr[$j] eq $v_arr[$j]) {
                        $self->found_match($match_flag, $j, \%coords);
                        $self->stats_incr('match', 1);
                        $match_flag = 1;
                    # mismatch
                    } else {
                        $self->stats_incr('mismatch', 1);
                        $match_flag = 0;
                    }
                }
            }
            $self->stats_incr('bp', scalar(@e_arr));
            $self->stats_incr('alignments', 1);
        }
        
        $i++;
    }
}

=head2 found_match

  Arg[4]      : Boolean $match_flag - flag indicating if last bp was a match
  Arg[5]      : Int $j - current bp position in the alignment
  Arg[7]      : Hashref $coords - alignment coordinates and strand from blastz
                output
  Description : Populates a datastructure describing blocks of alignment
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub found_match {
    my ($self, $match_flag, $j, $coords) = @_;

    my $id = $self->id;
    my $align = $self->get_stats('alignments');
    my $V_chr = $self->seq_region_name;

    # last position was a match
    if ($match_flag) {
        # adjust align block end
        if ($self->{'_match'}->{$V_chr}->{$id}->[$align]) {
            my $c = scalar(@{ $self->{'_match'}->{$V_chr}->{$id}->[$align] }) - 1;
            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[1] =
                $coords->{'e_start'} + $j - $self->get_stats('e_gap');
            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[3] =
                $coords->{'v_start'} + $j - $self->get_stats('v_gap');
        }
    
    # last position was a non-match
    } else {
        # start a new align block
        push @{ $self->{'_match'}->{$V_chr}->{$id}->[$align] }, [
            $coords->{'e_start'} + $j - $self->get_stats('e_gap'),
            $coords->{'e_start'} + $j - $self->get_stats('e_gap'),
            $coords->{'v_start'} + $j - $self->get_stats('v_gap'),
            $coords->{'v_start'} + $j - $self->get_stats('v_gap'),
            $coords->{'strand'},
            $coords->{'v_id'},
        ];
    }
}


=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub adjust_coords {
    my ($self, $e_start, $e_end, $V_coords) = @_;
    my $V_chr = $self->seq_region_name;
    my $id = $self->id;

    for (my $align = 0; $align < scalar(@{ $self->{'_match'}->{$V_chr}->{$id} || [] }); $align++) {
        for (my $c = 0; $c < scalar(@{ $self->{'_match'}->{$V_chr}->{$id}->[$align] || []}); $c++) {
            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[0] += $e_start - 1;
            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[1] += $e_start - 1;

            # forward strand match
            if ($self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[4] == 1) {
                $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[2] += $V_coords->{$self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[5]}->[0] - 1;
                $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[3] += $V_coords->{$self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[5]}->[0] - 1;
            
            # reverse strand match
            } else {
                my $tmp_start = $V_coords->{$self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[5]}->[1] - $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[3] + 1;

                $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[3] = $V_coords->{$self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[5]}->[1] - $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[2] + 1;

                $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[2] = $tmp_start;
            }

            # sanity check: aligned region pairs must have same length
            my $e_len = $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[1] - $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[0];
            my $v_len = $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[3] - $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[2];
            $self->support->log_warning("Length mismatch: $e_len <> $v_len in block $id, alignment $align, stretch $c\n", 2) unless ($e_len == $v_len);
        }
    }
}

=head2 filter_overlaps

  Description : Filters overlapping target (i.e. Vega) sequences in alignments.
                Longer alignments are preferred.
  Return type : none
  Exceptions  : none
  Caller      : internal

=cut

sub filter_overlaps {
    my $self = shift;
    my $id = $self->id;

    foreach my $V_chr (sort keys %{ $self->{'_match'} }) {
        # rearrange the datastructure so that we can find overlaps
        my $coord_check;
        foreach my $id (keys %{ $self->{'_match'}->{$V_chr} }) {
            for (my $align = 0; $align < scalar(@{ $self->{'_match'}->{$V_chr}->{$id} }); $align++) {
                for (my $c = 0; $c < scalar(@{ $self->{'_match'}->{$V_chr}->{$id}->[$align] || []}); $c++) {
                    push @{ $coord_check }, [
                        $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[0],
                        $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[1],
                        $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[2],
                        $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[3],
                        $id,
                        $align,
                        $c,
                    ];
                }
            }
        }
        
        my @e_sort = sort { $a->[0] <=> $b->[0] } @{ $coord_check };
        my @v_sort = sort { $a->[2] <=> $b->[2] } @{ $coord_check };

        # sanity check: Ensembl alignments must not overlap (axtBest should
        # guarantee that)
        my $last;
        foreach my $c (@e_sort) {
            $self->support->log_warning("Overlapping Ensembl alignment at ".join(':', $V_chr, $c->[0], $c->[1])." (last_end ".$last->[1].")\n", 1) if ($last and $c->[0] <= $last->[1]);
            $last = $c;
        }

        # now filter Vega overlaps
        my @seen;
        $last = undef;
        foreach my $c (@v_sort) {
            if ($last and $c->[2] <= $last->[3]) {
                $self->support->log_verbose("Overlapping Vega alignment at ".join(':', $V_chr, $c->[2], $c->[3])." (last_end ".$last->[3].")\n", 1);

                # if last alignment was longer, delete this one
                if ($last->[3]-$last->[2] > $c->[3]-$c->[2]) {
                    undef $self->{'_match'}->{$V_chr}->{$c->[4]}->[$c->[5]]->[$c->[6]];

                # if last alignment was shorter, trace back and delete all
                # overlapping shorter alignments
                } else {
                    foreach my $s (@seen) {
                        # earlier alignment still overlapping
                        if ($c->[2] <= $s->[3]) {
                            # earlier alignment shorter -> delete it
                            if ($s->[3]-$s->[2] < $c->[3]-$c->[2]) {
                                undef $self->{'_match'}->{$V_chr}->{$s->[4]}->[$s->[5]]->[$s->[6]];

                            # this alignment shorter -> delete it
                            } else {
                                undef $self->{'_match'}->{$V_chr}->{$c->[4]}->[$c->[5]]->[$c->[6]];
                                $last = $s;
                                last;
                            }
                        } else {
                            $last = $s;
                            last;
                        }
                    }                
                    
                    $last = $c;
                }
            }
            unshift @seen, $c;
            $last = $c unless ($last);
        }
    }
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub write_assembly {
    my ($self, $V_dba, $E_dbh, $E_sa) = @_;

    my $ensembl_chr_map = $self->support->get_ensembl_chr_mapping($V_dba, $self->support->param('assembly'));

    my $sth = $E_dbh->prepare(qq(
        INSERT INTO assembly (asm_seq_region_id, cmp_seq_region_id, asm_start,
            asm_end, cmp_start, cmp_end, ori)
        VALUES (?, ?, ?, ?, ?, ?, ?)
    ));

    $self->support->log("Adding assembly entries for alignments...\n");

    my $i;
    foreach my $V_chr (sort _by_chr_num keys %{ $self->{'_match'} }) {
        # get seq_region_id for Ensembl and Vega chromosome
        my $E_chr = $ensembl_chr_map->{$V_chr};
        my $V_sid = $E_sa->get_seq_region_id($E_sa->fetch_by_region('chromosome', $V_chr, undef, undef, undef, $self->support->param('assembly')));
        my $E_sid = $E_sa->get_seq_region_id($E_sa->fetch_by_region('chromosome', $E_chr, undef, undef, undef, $self->support->param('ensemblassembly')));

        foreach my $id (sort { $a <=> $b } keys %{ $self->{'_match'}->{$V_chr} }) {
            for (my $align = 0; $align < scalar(@{ $self->{'_match'}->{$V_chr}->{$id} }); $align++) {
                for (my $c = 0; $c < scalar(@{ $self->{'_match'}->{$V_chr}->{$id}->[$align] }); $c++) {
                    if ($self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]) {
                        $sth->execute(
                            $V_sid,
                            $E_sid,
                            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[2],
                            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[3],
                            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[0],
                            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[1],
                            $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[4],
                        );
                        $i++;
                    }
                }
            }
        }
    }
    $self->support->log("Done inserting $i entries.\n");
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub stats_incr {
    my ($self, $code, $incr) = @_;
    $self->{'_stats'}->{$code} += $incr;
    return $self->{'_stats'}->{$code};
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub init_stats {
    my ($self, @codes) = @_;
    foreach my $code (@codes) {
        $self->{'_stats'}->{$code} = 0;
    }
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub get_stats {
    my ($self, $code) = @_;
    return $self->{'_stats'}->{$code};
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub log_block_stats {
    my ($self, $indent) = @_;

    $self->support->log("Blastz alignment stats:\n", $indent);
    $self->support->log(sprintf(FMT2, "Alignments:", $self->get_stats('alignments')), $indent+1);
    if ($self->get_stats('alignments')) {
        $self->support->log(sprintf(FMT1,
            "Matches:",
            $self->get_stats('match'),
            $self->get_stats('match')/$self->get_stats('bp')*100),
        $indent+1);
        $self->support->log(sprintf(FMT1,
            "Mismatches:",
            $self->get_stats('mismatch'),
            $self->get_stats('mismatch')/$self->get_stats('bp')*100),
        $indent+1);
        $self->support->log(sprintf(FMT1,
            "Gaps:",
            $self->get_stats('gap'),
            $self->get_stats('gap')/$self->get_stats('bp')*100),
        $indent+1);
    }
    map { $self->stats_incr($_.'_total', $self->get_stats($_)) }
        qw(match mismatch gap bp);
}

=head2 

  Arg[1]      : 
  Example     : 
  Description : 
  Return type : 
  Exceptions  : 
  Caller      : 

=cut

sub log_overall_stats {
    my $self = shift;
    
    # blastz
    $self->support->log("\nOverall blastz alignment stats:\n");
    $self->support->log(sprintf(FMT1,
        "Matches:",
        $self->get_stats('match_total'),
        $self->get_stats('match_total')/$self->get_stats('bp_total')*100),
    1);
    $self->support->log(sprintf(FMT1,
        "Mismatches:",
        $self->get_stats('mismatch_total'),
        $self->get_stats('mismatch_total')/$self->get_stats('bp_total')*100),
    1);
    $self->support->log(sprintf(FMT1,
        "Gaps:",
        $self->get_stats('gap_total'),
        $self->get_stats('gap_total')/$self->get_stats('bp_total')*100),
    1);

    # alignments to be written to assembly table
    $self->support->log_verbose("\nAlignments that will be written to assembly table:\n");
    $self->support->log_verbose(sprintf(FMT3,
        qw(CHR BLOCK ALIGNMENT E_START E_END V_START V_END)),
    1);
    $self->support->log_verbose(('-'x63)."\n", 1);
    foreach my $V_chr (sort _by_chr_num keys %{ $self->{'_match'} }) {
        foreach my $id (sort { $a <=> $b } keys %{ $self->{'_match'}->{$V_chr} }) {
            for (my $align = 0; $align < scalar(@{ $self->{'_match'}->{$V_chr}->{$id} }); $align++) {
                for (my $c = 0; $c < scalar(@{ $self->{'_match'}->{$V_chr}->{$id}->[$align] }); $c++) {
                    if ($self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]) {
                        $self->support->log_verbose(sprintf(FMT4,
                            $V_chr,
                            $id,
                            $align+1,
                            @{ $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c] }),
                        1);
                        my $l = $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[1] - $self->{'_match'}->{$V_chr}->{$id}->[$align]->[$c]->[0];
                        $self->stats_incr('alignments_total', 1);
                        $self->stats_incr('short1_10_total', 1) if ($l < 11);
                        $self->stats_incr('short11_100_total', 1) if ($l > 10 and $l < 101);
                    }
                }
            }
        }
    }
    $self->support->log("\nAssembly entry stats:\n");
    $self->support->log(sprintf(FMT2,
        "Total alignment blocks:",
        $self->get_stats('alignments_total')),
    1);
    $self->support->log(sprintf(FMT2,
        "Alignments 1-10 bp:",
        $self->get_stats('short1_10_total')),
    1);
    $self->support->log(sprintf(FMT2,
        "Alignments 11-100 bp:",
        $self->get_stats('short11_100_total')),
    1);
}

=head2 _by_chr_num

  Example     : my @sorted = sort _by_chr_num qw(X, 6-COX, 14, 7);
  Description : Subroutine to use in sort for sorting chromosomes. Sorts
                numerically, then alphabetically
  Return type : values to be used by sort
  Exceptions  : none
  Caller      : internal

=cut

sub _by_chr_num {
    my @awords = split /-/, $a;
    my @bwords = split /-/, $b;

    my $anum = $awords[0];
    my $bnum = $bwords[0];

    if ($anum !~ /^[0-9]*$/) {
        if ($bnum !~ /^[0-9]*$/) {
            return $anum cmp $bnum;
        } else {
            return 1;
        }
    }
    if ($bnum !~ /^[0-9]*$/) {
        return -1;
    }

    if ($anum <=> $bnum) {
        return $anum <=> $bnum;
    } else {
        if ($#awords == 0) {
            return -1;
        } elsif ($#bwords == 0) {
            return 1;
        } else {
            return $awords[1] cmp $bwords[1];
        }
    }
}

=head2 AUTOLOAD

  Arg[1]      : (optional) String/Object - attribute to set
  Example     : # setting a attribute
                $self->attr($val);
                # getting the attribute
                $self->attr;
                # undefining an attribute
                $self->attr(undef);
  Description : lazy function generator for getters/setters
  Return type : String/Object
  Exceptions  : none
  Caller      : general

=cut

sub AUTOLOAD {
    my $self = shift;
    my $attr = our $AUTOLOAD;
    $attr =~ s/.*:://;
    return unless $attr =~ /[^A-Z]/;
    no strict 'refs';
    *{$AUTOLOAD} = sub {
        $_[0]->{'_data'}->{$attr} = $_[1] if (@_ > 1);
        return $_[0]->{'_data'}->{$attr};
    };
    $self->{'_data'}->{$attr} = shift if (@_);
    return $self->{'_data'}->{$attr};
}

1;

