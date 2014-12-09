package EnsemblAPI;

###################################
# Package to connect to Ensembl API
###################################


$VERSION = v0.0.1;

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;


sub getOrtholog {

#### Variables
# Connect API
my $reg			=	"Bio::EnsEMBL::Registry";
$reg->load_registry_from_db(
#   -host    	=> 'genobdd',
#   -user    	=> 'AUTOGRAPH2',
#   -pass    	=> 'Czbe8ULTtNUqK99Y',
#   -dbname  	=> 'ensembl_compara'
	-host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);
# To avoid losing conneection
$reg->set_reconnect_when_lost();


	my $id                  = shift;	    # ENSMBL ID
	my $targetspecies       = shift;	    #	<Target_species> Homo_sapiens || Canis_familiaris ||  Mus_musculus 
	my $orthtype            = "one2one";
	my $verbosity           = shift;
	$verbosity       		= 0 unless defined ($verbosity);


    # Input file
    if ($verbosity>0){ print STDERR "Looking for $id ortholog in $targetspecies...\n";}


	my $member_adaptor = $reg->get_adaptor('Multi', 'compara', 'Member');
	#print Dumper $member_adaptor;
	
	# Get member
	print STDERR "in sub: $id\n";
	my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$id);

	# Test member exist
	if (!defined  $member){
		return "NA";
	}

	# then you get the homologies where the member is involved
	my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');
	
	# get homologies with target species
	# That will return a reference to an array with all homologies (orthologues in
	# other species and paralogues in the same one)
	# Then for each homology, you can get all the Members implicated
	my $homologies = $homology_adaptor->fetch_all_by_Member_paired_species($member,$targetspecies);
 	
 	# Check if there are any elements in the referenced array:
	if (! @$homologies) {
		return "NA";
	} 
	else {
		foreach my $homology (@{ $homologies }){

			my $description = $homology->description;
			$description =~ s/possible_ortholog/one2many/;
			if ($description =~ /$orthtype/){
				foreach my $member (@{$homology->get_all_Members}){

# 					print Dumper $member;
					# member ortholog is done at the protein level
					# example : ENSCAFP00000002515.....Transcript:ENSCAFT00000002708 Gene:ENSCAFG00000001726 Chr:14 Start:8728365 End:8735331....
					if ($member->description !~ /$id/){				
						my @tmp	=	split(' ',$member->description);
						$tmp[1] =~	s/Gene://;
						#print "$id	$tmp[1]	$description\n";
						return $tmp[1];
					}
				}			
			}
		}
	}	
}

1;
