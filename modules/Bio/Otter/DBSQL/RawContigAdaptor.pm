package Bio::Otter::DBSQL::RawContigAdaptor;

use strict;

use Bio::EnsEMBL::DBSQL::RawContigAdaptor;

use vars qw(@ISA);

@ISA = qw ( Bio::EnsEMBL::DBSQL::RawContigAdaptor);


# Need to override fetch_filled_by_dbIDs in Bio::EnsEMBL::DBSQL::RawContigAdaptor
# because it makes Clones without going throught the CloneAdaptor

sub fetch_filled_by_dbIDs {
   my ($self,@ids) = @_;

   my $result = $self->SUPER::fetch_filled_by_dbIDs(@ids);

   foreach my $key (keys %$result) {
      my $clone = $result->{$key}->clone;

      $clone = $self->db->get_CloneAdaptor->annotate_clone($clone);  

   }

   return $result;
}


1;

	





