package Bio::Otter::DBSQL::DBAdaptor;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::Otter::DBSQL::AuthorGroupAdaptor;
@ISA = qw(Bio::EnsEMBL::DBSQL::DBAdaptor);

=head2 new

  Description: see Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($dataset)  = $self->rearrange([qw(DATASET)],@args);
  $self->dataset($dataset) if $dataset;
  return $self;
}

sub dataset {
  my ($self,$dataset) = @_;
  if (defined($dataset)) {
    $self->{_dataset} = $dataset;
  }
  return  $self->{_dataset};
}

sub begin_work {
    my $self = shift;
    $self->db_handle->begin_work();
}

sub commit {
    my $self = shift;
    $self->db_handle->commit();
}

sub rollback {
    my $self = shift;
    $self->db_handle->rollback();
}

sub get_AuthorGroupAdaptor {
  my $self = shift;
  if ( !exists $self->{'AuthorGroup'} ){
	 my $ad=Bio::Otter::DBSQL::AuthorGroupAdaptor->new($self);
	 $self->{'AuthorGroup'}=$ad;
  }
  return $self->{'AuthorGroup'};
}

=head2 get_available_adaptors

  Example     : my %object_adaptors = %{ $dbadaptor->get_available_adaptors };
  Description : returns a lookup hash of object adaptors for this DBAdaptor
  Return type : Hashref
  Exceptions  : none
  Caller      : Bio::EnsEMBL::Utils::ConfigRegistry

=cut

sub get_available_adaptors{
  my %pairs =  ( 
    'Author'               => 'Bio::Otter::DBSQL::AuthorAdaptor',
    'Clone'                => 'Bio::Otter::DBSQL::AnnotatedCloneAdaptor',
    'CloneInfo'            => 'Bio::Otter::DBSQL::CloneInfoAdaptor',
    'CloneLock'            => 'Bio::Otter::DBSQL::CloneLockAdaptor',
    'CloneRemark'          => 'Bio::Otter::DBSQL::CloneRemarkAdaptor',
    'Evidence'             => 'Bio::Otter::DBSQL::EvidenceAdaptor',
    'Gene'                 => 'Bio::Otter::DBSQL::AnnotatedGeneAdaptor',
    'GeneInfo'             => 'Bio::Otter::DBSQL::GeneInfoAdaptor',
    'CurrentGeneInfo'      => 'Bio::Otter::DBSQL::CurrentGeneInfoAdaptor',
    'GeneName'             => 'Bio::Otter::DBSQL::GeneNameAdaptor',
    'GeneSynonym'          => 'Bio::Otter::DBSQL::GeneSynonymAdaptor',
    'GeneRemark'           => 'Bio::Otter::DBSQL::GeneRemarkAdaptor',
    'Keyword'              => 'Bio::Otter::DBSQL::KeywordAdaptor',
    'MetaContainer'        => 'Bio::Otter::DBSQL::MetaContainer',
    'StableId'             => 'Bio::Otter::DBSQL::StableIdAdaptor',
    'Transcript'           => 'Bio::Otter::DBSQL::AnnotatedTranscriptAdaptor',
    'TranscriptClass'      => 'Bio::Otter::DBSQL::TranscriptClassAdaptor',
    'TranscriptInfo'       => 'Bio::Otter::DBSQL::TranscriptInfoAdaptor',
    'CurrentTranscriptInfo'=> 'Bio::Otter::DBSQL::CurrentTranscriptInfoAdaptor',
    'TranscriptRemark'     => 'Bio::Otter::DBSQL::TranscriptRemarkAdaptor',

    'Analysis'             => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
    'ArchiveStableId'      => 'Bio::EnsEMBL::DBSQL::ArchiveStableIdAdaptor',
    'Attribute'            => 'Bio::EnsEMBL::DBSQL::AttributeAdaptor',
    'AssemblyExceptionFeature' => 'Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor',
    'AssemblyMapper'       => 'Bio::EnsEMBL::DBSQL::AssemblyMapperAdaptor',
    # 'Blast'                => 'Bio::EnsEMBL::External::BlastAdaptor',
    # 'MetaContainer'        => 'Bio::EnsEMBL::DBSQL::MetaContainer',
    'CoordSystem'          => 'Bio::EnsEMBL::DBSQL::CoordSystemAdaptor',
    'CompressedSequence'   => 'Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor',
    'DBEntry'              => 'Bio::EnsEMBL::DBSQL::DBEntryAdaptor',
    'DnaAlignFeature'      => 'Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor',
    'DensityFeature'       => 'Bio::EnsEMBL::DBSQL::DensityFeatureAdaptor',
    'DensityType'          => 'Bio::EnsEMBL::DBSQL::DensityTypeAdaptor',
    'Exon'                 => 'Bio::EnsEMBL::DBSQL::ExonAdaptor',
    # 'Gene'                 => 'Bio::EnsEMBL::DBSQL::GeneAdaptor',
    'KaryotypeBand'        => 'Bio::EnsEMBL::DBSQL::KaryotypeBandAdaptor',
    'Marker'               => 'Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor',
    'MarkerFeature'        => 'Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor',
    'MetaCoordContainer'   => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
    'MiscSet'              => 'Bio::EnsEMBL::DBSQL::MiscSetAdaptor',
    'MiscFeature'          => 'Bio::EnsEMBL::DBSQL::MiscFeatureAdaptor',
    'PredictionTranscript' => 'Bio::EnsEMBL::DBSQL::PredictionTranscriptAdaptor',
    'PredictionExon'       => 'Bio::EnsEMBL::DBSQL::PredictionExonAdaptor',
    'ProteinFeature'       => 'Bio::EnsEMBL::DBSQL::ProteinFeatureAdaptor',
    'ProteinAlignFeature'  => 'Bio::EnsEMBL::DBSQL::ProteinAlignFeatureAdaptor',
    'SNP'                  => 'Bio::EnsEMBL::DBSQL::ProxySNPAdaptor',
    'QtlFeature'           => 'Bio::EnsEMBL::Map::DBSQL::QtlFeatureAdaptor',
    'Qtl'                  => 'Bio::EnsEMBL::Map::DBSQL::QtlAdaptor',
    'RepeatConsensus'      => 'Bio::EnsEMBL::DBSQL::RepeatConsensusAdaptor',
    'RepeatFeature'        => 'Bio::EnsEMBL::DBSQL::RepeatFeatureAdaptor',
    'Sequence'             => 'Bio::EnsEMBL::DBSQL::SequenceAdaptor',
    'SimpleFeature'        => 'Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor',
    'Slice'                => 'Bio::EnsEMBL::DBSQL::SliceAdaptor',
    'SupportingFeature'    => 'Bio::EnsEMBL::DBSQL::SupportingFeatureAdaptor',
    'TranscriptSupportingFeature'    => 'Bio::EnsEMBL::DBSQL::TranscriptSupportingFeatureAdaptor',
    # 'Transcript'           => 'Bio::EnsEMBL::DBSQL::TranscriptAdaptor',
    'Translation'          => 'Bio::EnsEMBL::DBSQL::TranslationAdaptor'
  );
  return (\%pairs);
}

sub get_MetaContainer {
  my $self = shift;
  $self->get_MetaContainerAdaptor();
}



1;

