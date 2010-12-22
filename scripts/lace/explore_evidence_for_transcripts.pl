#!/usr/bin/env perl

use warnings;

use strict;
use Carp;
use IO::String;
use List::Util qw(min max);

use Bio::Otter::Lace::Defaults;
use Bio::Otter::Utils::MM;

use Bio::SeqIO;
use Bio::EnsEMBL::Pipeline::SeqFetcher;
use Bio::Factory::EMBOSS;       # EMBOSS needs to be on PATH - /software/pubseq/bin/EMBOSS-5.0.0/bin
                                # To verify, check that 'wossname water' runs successfully
use Bio::AlignIO;

use Hum::Pfetch;

my $opts;

{
    my $dataset_name = undef;
    $opts = {
        total => 0,
        quiet => 0,
        verbose => 0,
        dump_seq => 0,
        dump_features => 0,
        evi_type => undef,
        max_length => undef,
        max_features => undef,
        logic_names => undef,
    };

    my $usage = sub { exec('perldoc', $0) };
    # This do_getopt() call is needed to parse the otter config files
    # even if you aren't giving any arguments on the command line.
    Bio::Otter::Lace::Defaults::do_getopt(
        'h|help!'       => $usage,
        'dataset=s'     => \$dataset_name,
        'quiet!'        => \$opts->{quiet},
        'verbose!'      => \$opts->{verbose},
        'total!'        => \$opts->{total},
        'dumpseq!'      => \$opts->{dump_seq},
        'dumpfeatures!' => \$opts->{dump_features},
        'type=s'        => \$opts->{evi_type},
        'maxlength=i'   => \$opts->{max_length},
        'maxfeatures=i' => \$opts->{max_features},
        'logicnames=s'  => \$opts->{logic_names},
        ) or $usage->();

    $usage->() unless $dataset_name;

    if (my $et = $opts->{evi_type}) {
        unless (   $et eq 'ncRNA' 
                || $et eq 'EST'
                || $et eq 'Protein'
                || $et eq 'cDNA'
                || $et eq 'Genomic'
            ) {
            croak "type must be one of EST,ncRNA,Protein,cDNA,Genomic";
        }
    }

    if ($opts->{quiet} and not $opts->{total}) {
        carp "Using -quiet but not -total - no output will be produced!";
    }

    if ($opts->{logic_names}) {
        $opts->{logic_names} = [ split(',', $opts->{logic_names}) ];
    }

    {
        my %evitype2logic = (
            'EST' => ['Est2genome_human','Est2genome_human_raw'],
            );

        if ($opts->{evi_type} and not $opts->{logic_names}) {
            $opts->{logic_names} = $evitype2logic{$opts->{evi_type}};
            if ($opts->{logic_names}) {
                carp("Logic names set to '", join(',', @{$opts->{logic_names}}),
                     "' for evidence type '", $opts->{evi_type}, "'");
            } else {
                carp("No logic name translation for evidence type '", $opts->{evi_type}, "'");
            }
        }
    }

    # Client communicates with otter HTTP server
    my $cl = Bio::Otter::Lace::Defaults::make_Client();
    
    # DataSet interacts directly with an otter database
    my $ds = $cl->get_DataSet_by_name($dataset_name);
    my $otter_dba = $ds->get_cached_DBAdaptor;

    my $where = "";
    my @args = ();
    if ($opts->{evi_type}) {
        $where = "AND e.type = ?";
        push @args, $opts->{evi_type};
    }
    my $list_transcripts = $otter_dba->dbc->prepare(qq{
        SELECT DISTINCT
                t.transcript_id
        FROM
                transcript t
           JOIN gene g USING (gene_id)
           JOIN evidence e ON t.transcript_id = e.transcript_id
           JOIN seq_region_attrib sra ON t.seq_region_id = sra.seq_region_id
        WHERE
                t.is_current = 1
            AND g.source = 'havana'
            -- Make sure it is on a writeable seq_region
            AND sra.attrib_type_id = (SELECT attrib_type_id FROM attrib_type WHERE code = 'write_access')
            AND sra.value = 1
            $where
        ORDER BY t.transcript_id
    });
    $list_transcripts->execute(@args);

    my $count = 0;
    while (my ($tid) = $list_transcripts->fetchrow()) {
        ++$count;
        reportf("verbose:ML:0", "TID: %10d", $tid );
        process_transcript($otter_dba, $tid, $opts);
    }
    reportf("all:ML:0", "Total: %d", $count) if $opts->{total};

}

my ($seq_str, $seq_str_io, $seqio_out);

sub setup_io {
    $seq_str_io = IO::String->new(\$seq_str);
    $seqio_out = Bio::SeqIO->new(-format => 'Fasta',
                                 -fh     => $seq_str_io );
}

my $fetcher;

# Warning - Bio::EnsEMBL::Pipeline::SeqFetcher spawns a pfetch each time to do the work
#          
# Therefore no longer used
#
sub pfetch_ensembl_pipeline {
    my ( $id ) = @_;
    $fetcher ||= Bio::EnsEMBL::Pipeline::SeqFetcher->new;
    my $seq = $fetcher->run_pfetch($id);
    carp sprintf "Cannot pfetch '%s'!\n", $id unless $seq;
    return $seq;
}

sub pfetch {
    my ( $id ) = @_;
    my ($hum_seq) = Hum::Pfetch::get_Sequences($id);
    unless ($hum_seq) {
        carp sprintf "Cannot pfetch '%s'!\n", $id;
        return undef;
    }
    my $seq = Bio::Seq->new(
        -seq        => $hum_seq->sequence_string,
        -display_id => $hum_seq->name,
        );
    return $seq;
}

my $mm;

sub get_accession_type {
    my $name = shift;

    $mm ||= Bio::Otter::Utils::MM->new;

    my $accession_types = $mm->get_accession_types([$name]);
    my $at = $accession_types->{$name};
    return @$at;
}

my ($factory, $comp_app);

sub get_comp_app {
    $factory ||= Bio::Factory::EMBOSS->new();
    $comp_app ||= $factory->program('water');
    return $comp_app;
}

sub process_align_features {
    my $features_ref = shift;
    my $opts = shift;

    return undef if $opts->{max_features} and scalar(@$features_ref) > $opts->{max_features};

    my %by_dbid = ();
    my %by_hseqname = ();
    my @ranked = ();
    my %fingerprints = ();
    my %rank_for_feature = ();
    my %by_hseqname_ranked = ();
    my %hseq_anl_acc = ();
    my @flat_hseq_anl_acc = ();

    my $rank = 0;
    my $count = 0;

    my @features = sort { $b->score      <=> $a->score
                                          ||
                          $b->percent_id <=> $a->percent_id } @{$features_ref};

  FEATURE: foreach my $feature ( @features ) {

        my $hseqname = $feature->hseqname;

        if ($opts->{evi_type} and not $opts->{logic_names}) {
            my ($f_type, $f_ver) = get_accession_type($hseqname);
            $f_type ||= "";
            if ($f_type ne $opts->{evi_type}) {
                reportf("normal:PF:1", "Skipping %s, type is '%s'", $hseqname, $f_type);
                next FEATURE;
            }
        }

        my $dbid = $feature->dbID;
        $by_dbid{$dbid} = $feature;

        my $anid = $feature->analysis->dbID;

        push @{$by_hseqname{$hseqname}}, $feature;

        # Should I look at ordering and positioning wrt exons?
        #
        my $acc = $hseq_anl_acc{$hseqname}->{$anid};
        unless (defined $acc) {
            $acc = $hseq_anl_acc{$hseqname}->{$anid} = { hseqname    => $hseqname,
                                                         analysis_id => $anid, 
                                                         dbid        => $dbid,
                                                         logic_name  => $feature->analysis->logic_name, 
                                                         count => 0, score => 0, length => 0, identical => 0 };
            push @flat_hseq_anl_acc, $acc;
        }
        ++$acc->{count};
        $acc->{score}  += $feature->score;
        $acc->{length} += $feature->length;
        $acc->{identical} += ($feature->length * $feature->percent_id / 100.0);
        $acc->{percent_id} = $acc->{identical} / $acc->{length} * 100.0;

        my $fp = make_align_fingerprint($feature);
        my $f_rank;
        if (defined $fingerprints{$fp}) {
            $f_rank = $fingerprints{$fp};
        } else {
            $f_rank = $fingerprints{$fp} = $rank;
            $rank++;
        }

        push @{$ranked[$f_rank]}, $feature;

        $rank_for_feature{$dbid} = $f_rank;

        reportf("all:PF:1", "%d %s %s", $dbid, $hseqname, $fp) if $opts->{dump_features};
        ++$count;
    }

    foreach my $i ( 0..$#ranked ) {
        map { push @{$by_hseqname_ranked{$_->hseqname}->{$i}}, $_ } @{$ranked[$i]};
    }

    my $by_total_hseqname_anl = [ sort { $b->{score} <=> $a->{score}
                                                      ||
                                         $b->{analysis_id} <=> $a->{analysis_id} }
                                  @flat_hseq_anl_acc ];

    reportf("verbose:PF:1", "Processed %d features with %d fingerprints", $count, $rank);

    return {
        _fingerprints       => \%fingerprints,
        _by_dbid            => \%by_dbid,
        _by_hseqname        => \%by_hseqname,
        _by_hseqname_ranked => \%by_hseqname_ranked,
        _ranked             => \@ranked,
        _rank_for_feature   => \%rank_for_feature,
        _acc_by_hseqname_anl=> \%hseq_anl_acc,
        _by_total_hseqname_anl => $by_total_hseqname_anl,
    };
}

sub make_align_fingerprint {
    my $feature = shift;
    return sprintf("%d-%d (%+d)\t=> %d-%d (%+d)\t: %d %.1f %s",
                   $feature->hstart(), $feature->hend(), $feature->hstrand(),
                   $feature->start(),  $feature->end(),  $feature->strand(),
                   $feature->score(), $feature->percent_id(), $feature->cigar_string()
        );
}

sub feature_augment {
    my ($f_list, $hseqname_seen, $n_list) = @_;
    my $n = 0;
    foreach my $f ( @$n_list ) {
        my $hseqname    = $f->hseqname;
        my $fingerprint = make_align_fingerprint($f);
        my $key         = join(':', $hseqname, $fingerprint);

        next if $hseqname_seen->{$key};

        push(@$f_list, $f);
        $hseqname_seen->{$key} = 1;
        ++$n;
    }
    return $n;
}

my ($transcript_adaptor, $gene_adaptor, $slice_adaptor, $dna_feature_adaptor);
my ($pipe_dba, $p_slice_adaptor);

sub process_transcript {
    my ($otter_dba, $tid, $opts) = @_;

    $transcript_adaptor ||= $otter_dba->get_TranscriptAdaptor;
    $gene_adaptor       ||= $otter_dba->get_GeneAdaptor;
    $slice_adaptor      ||= $otter_dba->get_SliceAdaptor;

    $pipe_dba        ||= Bio::Otter::Lace::PipelineDB::get_pipeline_DBAdaptor($otter_dba);
    $p_slice_adaptor ||= $pipe_dba->get_SliceAdaptor;

    my $td = $transcript_adaptor->fetch_by_dbID($tid);
    if ($td) {

        reportf("normal:PT:0", "Transcript %s, strand %d", $td->stable_id, $td->strand);

        #TEMP debug hook
        if ($td->stable_id eq 'OTTHUMT00000109016') {
            $DB::single = 1;
        }

        my $gene = $gene_adaptor->fetch_by_transcript_id($tid);

        # my $slice = $slice_adaptor->fetch_by_gene_stable_id($gene->stable_id); # not req'd?

        my $min_start = min($td->seq_region_start, $gene->seq_region_start);
        my $max_end   = max($td->seq_region_end,   $gene->seq_region_end);
        
        reportf("verbose:PT:1",
               "Fetching from pipeline DB: %s [%s] @ %s %s from %d to %d (g %d to %d) (ts %d to %d)",
               $gene->stable_id,
               $gene->external_name || '-',
               $gene->coord_system_name,
               $gene->seq_region_name,
               $min_start,
               $max_end,
               $gene->seq_region_start,
               $gene->seq_region_end,
               $td->seq_region_start,
               $td->seq_region_end,
            );

        my $p_slice = $p_slice_adaptor->fetch_by_region($gene->coord_system_name,
                                                        $gene->seq_region_name,
                                                        $gene->start,
                                                        $gene->end);

        # Is this interesting in practice?
        my $genes = $p_slice->get_all_Genes;
        foreach my $gene (@$genes) {
            reportf("verbose:PT:1", "Got pipeline gene %s [%s]", $gene->dbID, $gene->display_id);
        }
        
        my $p_features;
        if ($opts->{logic_names}) {
            $p_features = [];
            my $p_f_seen = {};
            foreach my $logic_name ( @{$opts->{logic_names}} ) {
                feature_augment($p_features, $p_f_seen, $p_slice->get_all_DnaAlignFeatures($logic_name));
            }
        } else {
            $p_features = $p_slice->get_all_DnaAlignFeatures();
        }
        report("verbose:PT:1", "Got ", scalar @$p_features, " features");

        my $dna_align_features = process_align_features($p_features, $opts);

        my $tf = $dna_align_features->{_by_total_hseqname_anl}->[0];
        if ($tf) {
            reportf("normal:PT:1", "Top feature: %s %d [Score:%d Length:%d %%ID:%d]",
                    $tf->{hseqname}, $tf->{dbid}, $tf->{score}, $tf->{length}, $tf->{percent_id});
        }
        
        if ($opts->{dump_seq}) {
            setup_io() unless $seqio_out;
            $seq_str_io->truncate(0);   # reset to start of $seq_str
            $seqio_out->write_seq($td->seq);
            print $seq_str;
        }

        my $all_dna_features;

        my @evidence = @{$td->evidence_list};
        EVIDENCE: foreach my $evi (@evidence) {

            if ($opts->{evi_type}) {
                next EVIDENCE unless $evi->type eq $opts->{evi_type};
            }

            my $e_name = $evi->name;
            my @at = get_accession_type($e_name);
            reportf("verbose:PT:1", "Evidence: %s - %s [%s %s]", $e_name, $evi->type, $at[0], $at[1]);
            if ($e_name ne $at[1] and $at[1] =~ m/^$e_name/) {
                reportf("verbose:PT:1", "Adding version, %s => %s", $e_name, $at[1]);
                $e_name = $at[1];
            }
            my ($prefix, $short_name) = $e_name =~ m/^(\w+)?:?([\w\.]+)$/ ;
            $short_name ||= $e_name;

            my $hit_feature = undef;
            if (my $flist = $dna_align_features->{_by_hseqname_ranked}->{$short_name}) {
                my $top_rank = (sort keys %$flist)[0];
                my $features = $flist->{$top_rank};
                $hit_feature = $features->[0];
                unless ($opts->{quiet}) {
                    my $n_features = scalar(@$features);
                    reportf("normal:PT:2", "Found feature match for %s at rank %d with %d features, logic %s",
                           $short_name, $top_rank, $n_features, $hit_feature->analysis->logic_name );
                    reportf("verbose:PT:2", "Feature fingerprint:  %s", make_align_fingerprint($hit_feature));

                    my $acc_feats = $dna_align_features->{_acc_by_hseqname_anl}->{$short_name};
                    if ($acc_feats) {
                        my $max_ann = (sort {$b <=> $a} keys %$acc_feats)[0];
                        my $af = $acc_feats->{$max_ann};
                        reportf("normal:PT:2", "Acc feature length: %d", $af->{length});
                    }
                }
            } else {
                report("normal:PT:2", "No feature match for $short_name");

                # Do a little more digging 

              FEATURE_HUNT: {

                  # First check it's not hiding under a different logic name, if we previously narrowed
                  if ($opts->{logic_names}) {
                      $all_dna_features ||= $p_slice->get_all_DnaAlignFeatures();
                      report("verbose:PT:2", "Got ", scalar @$all_dna_features, " features using all logic names");
                      my @match_features = grep { $_->hseqname eq $short_name } @$all_dna_features;
                      if (@match_features) {
                          reportf("normal:PT:2",
                                  "Found %d exact matches using all logic names", scalar(@match_features));
                          if ($opts->{verbose}) {
                              foreach my $feature (@match_features) {
                                  reportf("verbose:PT:3", "%s : %s",
                                          $feature->hseqname, $feature->analysis->logic_name);
                              }
                          }
                          last FEATURE_HUNT;
                      }
                      # This stanza shouldn't get hit now that we add version to $e_name above where necessary
                      @match_features = grep { $_->hseqname =~ /^$short_name/ } @$all_dna_features;
                      if (@match_features) {
                          reportf("normal:PT:2",
                                  "Found %d partial matches using all logic names", scalar(@match_features));
                          if ($opts->{verbose}) {
                              foreach my $feature (@match_features) {
                                  reportf("verbose:PT:3", "%s : %s",
                                          $feature->hseqname, $feature->analysis->logic_name);
                              }
                          }
                          last FEATURE_HUNT;
                      }
                  }

                  # Next lookup by name
                  $dna_feature_adaptor ||= $pipe_dba->get_DnaAlignFeatureAdaptor;
                  my $dna_features_by_name = $dna_feature_adaptor->fetch_all_by_hit_name($short_name);
                  report("verbose:PT:2", "Got ", scalar @$dna_features_by_name, " features by name");
                  if (@$dna_features_by_name) {
                      reportf("normal:PT:2", "Found %d matches by name in entire DB", scalar(@$dna_features_by_name));
                      last FEATURE_HUNT;
                  }

                  # This stanza shouldn't get hit now that we add version to $e_name above where necessary
                  $dna_features_by_name = $dna_feature_adaptor->fetch_all_by_hit_name_unversioned($short_name);
                  report("verbose:PT:2", "Got ", scalar @$dna_features_by_name, " features by name_unversioned");
                  if (@$dna_features_by_name) {
                      reportf("normal:PT:2",
                              "Found %d matches by name_unversioned in entire DB", scalar(@$dna_features_by_name));
                      last FEATURE_HUNT;
                  }

                  report("normal:PT:2", "No luck looking for feature :-(");

                } # FEATURE_HUNT
            }

            my $evi_seq = pfetch($e_name);
            next EVIDENCE unless $evi_seq;

            if (    $opts->{max_length}
                    and ($evi_seq->length > $opts->{max_length})
                    and ($td->seq->length > $opts->{max_length})
                ) {
                my $msg = sprintf("Transcript %s and evidence %s both too long, skipping",
                                  $td->stable_id, $e_name);
                report("all:PT:0", $msg);
                carp $msg;
                next EVIDENCE;
            }

            # Compare them
            my $comp_app = get_comp_app();
            my $comp_fh = File::Temp->new();
            my $comp_outfile = $comp_fh->filename;

            my $rev_str = $evi_seq->seq;
            Bio::EnsEMBL::Utils::Sequence::reverse_comp(\$rev_str);
            my $rev_evi_seq = Bio::Seq->new(
                -seq        => $rev_str,
                -display_id => $evi_seq->display_id . ".rev",
                );

            my ($msg, $exp);
            if ($hit_feature) {
                if ($hit_feature->hstrand == $td->strand) {
                    $msg = "Feature hit and transcript strands match    - expecting forward match";
                    $exp = 1;
                } else {
                    $msg = "Feature hit and transcript strands mismatch - expecting reverse match";
                    $exp = -1;
                }
            } else {
                $msg = "No hit feature - no match prediction";
            }
            report("normal:PT:2", $msg);

            $comp_app->run({-asequence => $td->seq,
                            -bsequence => [$evi_seq, $rev_evi_seq],
                            -outfile   => $comp_outfile,
                            -aformat   => 'srspair',
                           });

            my $alnin = Bio::AlignIO->new(-format => 'emboss',
                                          -fh     => $comp_fh);

            my @aln_results;
            while ( my $aln = $alnin->next_aln ) {
                push @aln_results, $aln;
                # process the alignment -- these will be Bio::SimpleAlign objects
                my $name = $aln->get_seq_by_pos(2)->display_id;
                if ($opts->{verbose}) {
                    reportf("verbose:PT:3", "feature: %-15s : score: %8.1f, length: %5d, ident: %5.1f%%",
                           $name,
                           $aln->score,
                           $aln->length,
                           $aln->percentage_identity,
                        );
                } elsif (not $opts->{quiet}) {
                    reportf("normal:PT:3", "%s,%s,%s,%.1f,%d,%.1f,%d,%d",
                           $td->stable_id,
                           $name,
                           $evi->type,
                           $aln->score,
                           $aln->length,
                           $aln->percentage_identity,
                           $td->seq->length,
                           $evi_seq->length,
                        );
                }
            }

            my $top_hit = (sort { $b->score <=> $a->score } @aln_results)[0];
            if ($top_hit == $aln_results[0]) {
                report("verbose:PT:2", "Top match on forward");
                if ($exp) {
                    if ($exp == 1) {
                        report("normal:PT:2", "Match on forward as expected");
                    } else {
                        report("all:PT:2", "MATCH ON FORWARD, EXPECTING REVERSE");
                    }
                }
            } else {
                report("verbose:PT:2", "Top match on reverse");
                if ($exp) {
                    if ($exp == -1) {
                        report("normal:PT:2", "Match on reverse as expected");
                    } else {
                        report("all:PT:2", "MATCH ON REVERSE, EXPECTING FORWARD");
                    }
                }
            }
            
        } # EVIDENCE

        printf "\n" unless $opts->{quiet};

    } else {
        carp "Cannot retrieve transcript with id %d from adaptor";
    }
}

# Rational output of results
#
#   opts:  quiet  neither verbose
# all       Y      Y       Y
# normal    N      Y       Y
# verbose   N      N       Y

sub report_prefix {
    my $cond = shift;
    my ($level, $context, $indent) = split(':', $cond);

    if ($opts->{quiet}) {
        return undef unless $level eq 'all';
    }
    unless ($opts->{verbose}) {
        return undef if $level eq 'verbose';
    }

    if ($opts->{context}) {
        my $n = 4*($indent+1);
        my $prefix_format = "%-${n}s";
        return sprintf($prefix_format, $context);
    } else {
        return $indent ? ' ' x (4*$indent) : '';
    }
}

sub report {
    my $cond = shift;
    my $prefix = report_prefix($cond);
    return unless defined $prefix;

    print $prefix, @_, "\n";
}

sub reportf {
    my $cond = shift;
    my $prefix = report_prefix($cond);
    return unless defined $prefix;

    my $format = shift;
    my $msg = sprintf($format, @_);
    print $prefix, $msg, "\n";
}

__END__

=head1 NAME - explore_evidence_for_transcripts.pl

=head1 SYNOPSIS

explore_evidence_for_transcripts.pl -dataset <DATASET NAME> [-type <EVIDENCE_TYPE>] [-quiet] [-total]

=head1 DESCRIPTION

Explore evidence matching a transcript.

=head1 AUTHOR

Michael Gray B<email> mg13@sanger.ac.uk

