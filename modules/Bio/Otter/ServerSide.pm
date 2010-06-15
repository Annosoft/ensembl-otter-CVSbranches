package Bio::Otter::ServerSide;

use strict;
use warnings;
use Exporter;

use Bio::Otter::ServerQuery;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Otter::DBSQL::DBAdaptor;
use Bio::Vega::DBSQL::DBAdaptor;
use Bio::Otter::Author;
use Bio::Otter::Version;
use Bio::Otter::Lace::SatelliteDB;
use Bio::Otter::LogFile;

our @ISA       = qw(Exporter);
our @EXPORT    = qw();
our @EXPORT_OK = qw(
                    &server_log
                    &set_nph
                    &send_response
                    &error_exit
                    &odba_to_sdba
                    &get_mapper_dba
                    &get_slice
                    &fetch_mapped_features
                    &get_Author_from_CGI
                    &get_DBAdaptor_from_CGI_species
                    );
our %EXPORT_TAGS = (all => [qw(
                               server_log
                               set_nph
                               send_response
                               error_exit
                               odba_to_sdba
                               get_mapper_dba
                               get_slice
                               fetch_mapped_features
                               get_Author_from_CGI
                               get_DBAdaptor_from_CGI_species 
                               )
                            ],
                    );

sub server_log {
    my $line = shift @_;
    my $csn = $ENV{CURRENT_SCRIPT_NAME} || $0;
    print STDERR "[$csn] $line\n";
}

sub set_nph{
    my ($sq) = @_;
    error_exit('', 'I need a CGI object') unless $sq && UNIVERSAL::isa($sq, 'CGI');
    if (defined($ENV{SERVER_SOFTWARE})
    && ( $ENV{SERVER_SOFTWARE} =~ /libwww-perl-daemon/)) {
        # server_log('NOTE: Setting nph to 1');
        $sq->nph(1);
    }
}

sub send_response{
    my ($sq, $response, $wrap) = @_;

    server_log('Sending the response =====================');
    print $sq->header('text/plain') if $sq && UNIVERSAL::isa($sq, 'CGI');

    if($wrap) {
        print qq`<?xml version="1.0" encoding="UTF-8"?>\n`;
        print qq`<otter schemaVersion="$SCHEMA_VERSION" xmlVersion="$XML_VERSION">\n`;
    }

    print $response;

    if($wrap) {
        print "</otter>\n";
    }
}

sub error_exit {
    my ($sq, $reason) = @_;

    chomp($reason);

    send_response($sq, " <response>\n    ERROR:\n$reason\n </response>", 1);
    server_log("ERROR: $reason\n");

    exit(1);
}

sub odba_to_sdba {
    my ($sq, $odba, $pipehead, $metakey) = @_;

    server_log("called with: ".join(' ', map { "$_=".$sq->param($_) } $sq->param ));

    $metakey ||= '';

        # It may well be true that the caller
        # is interested in features from otter_db itself.
        # (This is NOT the default behaviour,
        #  so he has to specify it by setting metakey='.')

    if($metakey eq '.') {
        server_log("Connecting to the otter_db itself");
        return $odba;      # and ignore $pipehead flag for now
    }

    my $kind = 'satellite DB';

    if(! $metakey) {
        $metakey = $pipehead
            ? 'pipeline_db_head'
            : 'pipeline_db';
        $kind = 'pipeline DB'
    }

    server_log("connecting to the ".($pipehead?'NEW':'OLD')." schema $kind using [$metakey] meta entry...");

    my $class = $pipehead
        ? 'Bio::EnsEMBL::DBSQL::DBAdaptor'
        : 'Bio::Otter::DBSQL::DBAdaptor';

    my ($sdba, $sdb_options) =
        Bio::Otter::Lace::SatelliteDB::_get_DBAdaptor_and_options(
            $odba,
            $metakey,
            $class,
        );

    error_exit($sq, "Could not create satellite_db for '$metakey' in otter database")
        unless ($sdba);

    error_exit($sq, "No connection parameters for '$metakey' in otter database")
        unless ($sdb_options && keys %$sdb_options);

    $sdba->assembly_type($odba->assembly_type()) unless $pipehead;

    server_log("... with parameters: ".join(', ', map { "$_=".$sdb_options->{$_} } keys %$sdb_options ));

    return $sdba;
}

sub get_mapper_dba {
    my ($sq, $odba, $sdba, $pipehead) = @_;

    my $metakey    = $sq->param('metakey')    || ''; # defaults to pipeline
    my $target_asm = $sq->param('target_asm') || undef;

    if(!$metakey) {
        server_log("Working with pipeline_db directly, no remapping is needed.");
        return;
    } elsif($metakey eq '.') {
        server_log("Working with otter_db directly, no remapping is needed.");
        return;
    }

    my ($sdb_asm) = @{ $sdba->get_MetaContainer()->list_value_by_key('assembly.default') };
    $sdb_asm ||= $target_asm || 'UNKNOWN';

    if($sdb_asm eq 'Otter') {
        server_log("No remapping is needed, the satellite_db is another Otter database.");
        return;
    }

    if($pipehead) { # a head version of ensembl_db (non-pipeline genes)

            # Currently we keep assembly equivalency information in the pipeline_db_head seq_region_attrib.
            # Once otter_db is converted into new schema, we can keep this information there.
        my $pdba = odba_to_sdba($sq, $odba, 1, ''); # ensures we get new pipeline
        my $pipe_slice = get_slice($sq, $pdba, 1);
        my %asm_is_equiv = map { ($_->value() => 1) } @{ $pipe_slice->get_all_Attributes('equiv_asm') };

        if($asm_is_equiv{$sdb_asm}) { # we can simply rename instead of mapping

            my $chr_name = $sq->param('name');
            server_log("This chr is equivalent to '$chr_name' in our reference '$sdb_asm' assembly");
            $sq->param('csver', $sdb_asm);
            return;

        } else { # assemblies are guaranteed to differ!

            my $mapper_metakey = "mapper_db.${sdb_asm}";

                # see if 'mapper_db' meta link is defined:
            my ($mapper_val) = @{$odba->get_MetaContainer()->list_value_by_key($mapper_metakey)};
            if($mapper_val) {
                my $mdba;

                    # we may keep the mapping information in new otter|pipeline db or elsewhere:
                if($mapper_val eq '=otter_head') {
                    $mdba = $odba;
                } elsif($mapper_val eq '=pipeline_head') {
                    $mdba = $pdba;
                } else {
                    $mdba = odba_to_sdba($sq, $odba, 1, $mapper_metakey);
                }

                return ($mdba, $sdb_asm);
            } else {
                server_log("No '$mapper_metakey' defined in meta table => cannot map between assemblies => exiting");
                send_response($sq, '', 1);
                exit(0);
            }
        }

    } else { # non-head version, cannot guarantee correct mapping

        server_log("No remapping is being done, you're responsible for doing it on the client side.");
        return;
    }
}

sub get_slice { # codebase-independent version for scripts
    my ($sq, $dba, $pipehead) = @_;

    my $slice;

    my $cs    = $sq->param('cs') || 'chromosome';
    my $name  = $sq->param('name');
    my $type  = $sq->param('type');
    my $start = $sq->param('start');
    my $end   = $sq->param('end');

    if($pipehead) {

        my $strand= $sq->param('strand');
        my $csver = $sq->param('csver');
        if(!$csver && ($cs eq 'chromosome')) {
            $csver = 'Otter';
        }

            # The following statement ensures
            # that we use 'assembly type' as the chromosome name
            # only for Otter chromosomes.
            # Vega chromosomes will have simple names.
        my $segment_attr = (($cs eq 'chromosome') && ($csver eq 'Otter'))
			? 'type'
		    : 'name';
        my $segment_name = $sq->param($segment_attr);

        error_exit($sq, "$cs '$segment_attr' attribute not set ") unless $segment_name;

        $slice =  $dba->get_SliceAdaptor()->fetch_by_region(
            $cs,
	        $segment_name,
            $start,
            $end,
            $strand,
            $csver,
        );

    } else {

        error_exit($sq, "$cs 'name' attribute not set") unless $name;

        if($cs eq 'chromosome') {
            $start ||= 1;

            eval {
                my $chr_obj = $dba->get_ChromosomeAdaptor()->fetch_by_chr_name($name);
                $end ||= $chr_obj->length();
            };
            if($@) {
                server_log("Could not get chromosome '$name', returning an empty list");
                send_response($sq, '', 1);
                exit(0); # <--- this forces all the scripts to exit normally
            }

            $slice = $dba->get_SliceAdaptor()->fetch_by_chr_start_end(
                $name,
                $start,
                $end,
            );

            if($slice and ! @{ $slice->get_tiling_path() } ) {
                server_log('Could not get a slice, probably not (yet) loaded into satellite db');
                send_response($sq, '', 1);
                exit(0); # <--- this forces all the scripts to exit normally
            }
        } elsif($cs eq 'contig') {
            eval {
                $slice = $dba->get_RawContigAdaptor()->fetch_by_name(
                    $name,
                );
            };
            if($@) {
                server_log("Could not get contig '$name', returning an empty list");
                send_response($sq, '', 1);
                exit(0); # <--- this forces all the scripts to exit normally
            }

        } else {
            error_exit($sq, "Other coordinate systems are not supported");
        }

    }

    if(not $slice) {
        server_log('Could not get a slice, probably not (yet) loaded into satellite db');
        send_response($sq, '', 1);
        exit(0); # <--- this forces all the scripts to exit normally
    }

    return $slice;
}

sub fetch_mapped_features {
    my ($sq, $pipehead, $feature_name, $call_parms) = @_;

    my $fetching_method = shift @$call_parms;

    my $cs       = $sq->param('cs')      || 'chromosome';
    my $csver    = $sq->param('csver')   || undef;
    my $metakey  = $sq->param('metakey') || ''; # defaults to pipeline

    my $odba = get_DBAdaptor_from_CGI_species($sq, $OTTER_SPECIES, $pipehead);
    my $sdba = odba_to_sdba($sq, $odba, $pipehead, $metakey);

    my ($mdba, $sdba_asm) = get_mapper_dba($sq, $odba, $sdba, $pipehead);

    my $features = [];

    if($mdba) {
        server_log("Proceeding with mapping code");

        my $original_slice_on_mapper = get_slice($sq, $mdba, $pipehead);
        my $proj_segments_on_mapper = $original_slice_on_mapper->project( $cs, $sdba_asm );

        my $sa_on_target = $sdba->get_SliceAdaptor();

        foreach my $segment (@$proj_segments_on_mapper) {
            my $projected_slice_on_mapper = $segment->to_Slice();

            my $target_slice_on_target = $sa_on_target->fetch_by_region(
                $projected_slice_on_mapper->coord_system()->name(),
                $projected_slice_on_mapper->seq_region_name(),
                $projected_slice_on_mapper->start(),
                $projected_slice_on_mapper->end(),
                $projected_slice_on_mapper->strand(),
                $projected_slice_on_mapper->coord_system()->version(),
            );

            my $target_fs_on_target_segment
                = $target_slice_on_target->$fetching_method(@$call_parms);

            server_log('***** : '.scalar(@$target_fs_on_target_segment)." ${feature_name}s found on the slice");

            foreach my $target_feature (@$target_fs_on_target_segment) {

                if($target_feature->can('propagate_slice')) {
                    $target_feature->propagate_slice($projected_slice_on_mapper);
                } else {
                    $target_feature->slice($projected_slice_on_mapper);
                }

                if( my $transferred = $target_feature->transfer($original_slice_on_mapper) ) {
                    push @$features, $transferred;
                } else {
                    my $fname = sprintf( "%s [%d..%d]", 
                                        $target_feature->display_id(),
                                        $target_feature->start(),
                                        $target_feature->end() );
                    server_log("Could not transfer $feature_name $fname onto {$cs:$csver}");
                }
            }
        }

    } else {
        server_log("No mapping is needed, just fetching");

        my $original_slice = get_slice($sq, $sdba, $pipehead);

        $features = $original_slice->$fetching_method(@$call_parms);
    }

    return $features;
}

sub get_Author_from_CGI{
  my ($sq,$pipehead) = @_;
  error_exit('', 'I need a CGI object') unless $sq && UNIVERSAL::isa($sq, 'CGI');
  my $auth_name = $sq->param('author') || error_exit($sq, "Need author for this script...");
  my $email     = $sq->param('email')  || error_exit($sq, "Need email for this script...");
  my $author;
  unless($pipehead) {
	 $pipehead=0;
  }

  if ($pipehead) {
	 $author=Bio::Vega::Author->new(-name  => $auth_name,
											  -email => $email);
  }
  else {
	 $author    = Bio::Otter::Author->new(-name  => $auth_name,
													  -email => $email);
  }
  return $author;
}

sub get_DBAdaptor_from_CGI_species{
    my ($sq, $SPECIES, $pipehead) = @_;

    error_exit('', 'I need two arguments') unless $sq && $SPECIES;
    error_exit('', 'I need a CGI object') unless UNIVERSAL::isa($sq, 'CGI');

    ####################################################################
    # Check the dataset has been entered
    my $dataset = $sq->param('dataset') || error_exit($sq, "No dataset type entered.");

    # get the overriding dataset options from species.dat 
    my $dbinfo   = $SPECIES->{$dataset} || error_exit($sq, "Unknown data set $dataset");

    # get the defaults from species.dat
    my $defaults = $SPECIES->{'defaults'};

    ########## CODEBASE tricks ########################################
    my $headcode  = $dbinfo->{HEADCODE} || $defaults->{HEADCODE};
    $pipehead ||= $headcode;

    my $type     = $sq->param('type') || $dbinfo->{TYPE} || $defaults->{TYPE};

    ########## AND DB CONNECTION #######################################

    my $dbhost    = $dbinfo->{HOST}     || $defaults->{HOST};
    my $dbuser    = $dbinfo->{USER}     || $defaults->{USER};
    my $dbpass    = $dbinfo->{PASS}     || $defaults->{PASS};
    my $dbport    = $dbinfo->{PORT}     || $defaults->{PORT};
    my $dbname    = $dbinfo->{DBNAME}   ||
		error_exit($sq, "Failed opening otter database [No database name]");

    my $dnahost    = $dbinfo->{DNA_HOST}    || $defaults->{DNA_HOST};
    my $dnauser    = $dbinfo->{DNA_USER}    || $defaults->{DNA_USER};
    my $dnapass    = $dbinfo->{DNA_PASS}    || $defaults->{DNA_PASS};
    my $dnaport    = $dbinfo->{DNA_PORT}    || $defaults->{DNA_PORT};
    my $dna_dbname = $dbinfo->{DNA_DBNAME};
  

    my $adaptor_class = $pipehead
        ? ( $headcode
                ? 'Bio::Vega::DBSQL::DBAdaptor'     # headcode anyway, get the best adaptor
                : 'Bio::EnsEMBL::DBSQL::DBAdaptor'  # new pipeline of the old otter, get the minimal adaptor
          )
        : ( $headcode
                ? 'Bio::EnsEMBL::DBSQL::DBAdaptor'  # old pipeline of the new otter, get the minimal adaptor
                : 'Bio::Otter::DBSQL::DBAdaptor'    # oldcode anyway, get the best adaptor
        );

    my( $odba, $dnadb );

    server_log("OtterDB='$dbname' host='$dbhost' user='$dbuser' pass='$dbpass' port='$dbport'");
    eval {
       $odba = $adaptor_class->new( -host   => $dbhost,
                                    -user   => $dbuser,
                                    -pass   => $dbpass,
                                    -port   => $dbport,
                                    -dbname => $dbname);
    };
    error_exit($sq, "Failed opening otter database [$@]") if $@;
    if ($dna_dbname) {
        eval {
            $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dnahost,
                                                        -user   => $dnauser,
                                                        -pass   => $dnapass,
                                                        -port   => $dnaport,
                                                        -dbname => $dna_dbname);
        };
        error_exit($sq, "Failed opening dna database [$@]") if $@;
        $odba->dnadb($dnadb);
        
        server_log("Connected to dna database");
    }
    if(!$pipehead) {
        server_log("Assembly_type='" . $odba->assembly_type($type)."'");
    }
    return $odba;
}

1;

__END__

=head1 AUTHOR

Ana Code B<email> anacode@sanger.ac.uk
