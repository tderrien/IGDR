package Interactions;

=head1 DESCRIPTION

Interactions is defined in a pairwise manner between 2 features from 2 hashes
Typically, the case for lncRNA wrt to mRNA features
An interaction is defined by
- idA			= id of partnerA (lncRNA)
- idB			= id of partnerB (mRNA)
- type			= 	* 1 genic 
					* 0 intergenic
- subtype		= 	* exonic or intronic or contain
					* div(ergent) or conv(ergent) or same(strand) or unk(nown)
- orient		= 	sense or antisense or NA
- dist	 		= 	* 0 for overlapping feature
					* Integer corresponding to the dist between start and end of the 2 features
- oversize		= size of overlap in nt.
- transcript_biotypeA		= default "."
- transcript_biotypeB		= default "."

$partner{'id'}			= $tr2;
$partner{'size'}		= $tr2_size;
$partner{'source'}		= $h2{$tr2}->{"$source"};
$partner{'type'}		= 0 or 1 ;					
$partner{'dist'}		= Utils::getdistance($b1,$e1,$b2,$e2);					
$partner{'subtype'}		= intergenic_Subtype($s1, $s2, $b1, $b2);
$partner{'oversize'}	= 0;
			
=cut
# 
$VERSION = v0.0.1;

use warnings;
use strict;
use Utils;
use ExtractFromFeature;

use Scalar::Util;
use Data::Dumper;

$| = 1;


=head2 InteractionCollection

 Title   : InteractionCollection
 Usage   : my %interactions = InteractionCollection($h1, $h2, $window, $verbosity)
 Function: Returns the interactions between 2 hashes wrt to a window (in b.p) around the %h1
 Returns : hash
 Args    : 


=cut
sub InteractionCollection{


	my ($refh1, $refh2, $window, $filterout, $strandedmode, $verbosity)	= @_;
	my %h1			=	%{$refh1}; # parsing a hash in sub need dereference in shift : lncRNA
	my %h2			=	%{$refh2}; # parsing a hash in sub need dereference in shift : mRNA
	$window 		||= 1000; # default window around partnerA (lncRNA)
	$strandedmode 	||= 0; # 0 keep all interactions, 1 keep all interactions with same strand (sense) , -1 keep all interactions with different strand (sense)
	$filterout		||= "none"; # remove genic or intergenic interactions [ default : none (keep every interactions)]
	$verbosity 		||= 0;

	die "InteractionCollection: Filter '$filterout' not known\n" if ($filterout ne "none" && $filterout ne "genic" && $filterout ne "intergenic");
	
	# stranded mode 
	my $strandedmode = 0 ;
	
	my $nb_tx1 = keys( %h1 );
	my $i=0;
	
	my %h_inter = ();
	
	# Sort both hash by chr and transcript start ie startt
	foreach my $tr1 (sort {
        $h1{$a}->{'startt'} <=> $h1{$b}->{'startt'} }
               keys %h1) {

		# verbose on file1
	    if ($verbosity > 0){ Utils::showProgress( $nb_tx1, $i++, "Interaction Collection: ");}

		my $b1	= $h1{$tr1}->{"startt"}; my $e1	= $h1{$tr1}->{"endt"};	my $s1	= $h1{$tr1}->{"strand"};
		my $tr1_size	= ExtractFromHash::cumulSize($h1{$tr1}->{"feature"});

		$h_inter{$tr1}->{"chr"}			=   $h1{$tr1}->{"chr"};
		$h_inter{$tr1}->{"source"}		=   $h1{$tr1}->{"source"};
		$h_inter{$tr1}->{"startt"}		=   $b1;
		$h_inter{$tr1}->{"endt"}		=   $e1;
		$h_inter{$tr1}->{"strand"}      =   $s1;
		$h_inter{$tr1}->{"size"}		=   $tr1_size;
	

	    # Sort both hash by chr and transcript start ie startt
		foreach my $tr2 (sort { 
            $h2{$a}->{'startt'} <=> $h2{$b}->{'startt'} }
            keys %h2) {				              
            

            # define important variables for overlap
            my $b2	= $h2{$tr2}->{"startt"}; my $e2	= $h2{$tr2}->{"endt"}; my $s2	= $h2{$tr2}->{"strand"};

# 			print "$tr1 $b1 $e1 $s1 -- $tr2 $b2 $e2 $s2\n";
            
            # Go to the end if tx (B + window) > A 
            last if ( ($e1 + $window) < $b2 );
            
            # Next if  (B + window) < A 
        	next if ( ($b1 - $window) > $e2 );

			###													  ###
			############ We are in the overlap window ###############
			###													  ###
			
			# Hash containing info on the interaction
			my %partner = ();				
				
			# Compute infos for tx2
    	    my $tr2_size	= ExtractFromHash::cumulSize($h2{$tr2}->{"feature"});

 			$partner{'id'}		= $tr2;
 			$partner{'size'}	= $tr2_size;
			$partner{'transcript_biotype'}	= $h2{$tr2}->{"transcript_biotype"};
			   	    
			
			# Interaction type : genic opr intergnic
			my $isgenic = interaction_Type( $b1,$e1,$b2,$e2,$s1,$s2, $strandedmode);

			# Filterout 			
			next if ($filterout eq "genic" && $isgenic);
			next if ($filterout eq "intergenic" && !$isgenic);

    	    # if genic    	    
    	    if ($isgenic){
    	    	    	    	
 				$partner{'type'}		= 1 ;					
 				$partner{'dist'}		= 0;	
				$partner{'orient'}	= orientation($s1,$s2);
 								

	   	    	# **** EXON ***
				# Compute the overlap size for Exon
				my $cumul_overlapExon_size	=	ExtractFromFeature::intersectFeatures($h1{$tr1}->{'feature'}, $h2{$tr2}->{'feature'}, $s1, $s2, $strandedmode, $verbosity);
				# Compute the overlap fraction	Exon		
				my $fractionExonOver1	=	$cumul_overlapExon_size/$tr1_size;
				my $fractionExonOver2	=	$cumul_overlapExon_size/$tr2_size;
				
				if ($cumul_overlapExon_size > 0){
					
					$partner{'subtype'}		= "exonic";
					$partner{'oversize'}	= $cumul_overlapExon_size;					
				}
# 				printf ("Exonic :\t%s\t%s\t%d\t%.2f\t%.2f\n", $tr1,$tr2, $cumul_overlapExon_size, $fractionExonOver1, $fractionExonOver2);
			

    	    	# **** INTRON ***
				# Compute the overlap size for Intron
				my $cumul_overlapIntron_size	=	ExtractFromFeature::intersectExonsIntrons($h1{$tr1}->{'feature'}, $h2{$tr2}->{'feature'}, $s1, $s2, $strandedmode, $verbosity);
				# Compute the overlap fraction	Intron		
				my $fractionIntronOver1	=	$cumul_overlapIntron_size/$tr1_size;
				my $fractionIntronOver2	=	$cumul_overlapIntron_size/$tr2_size;	

				if ($cumul_overlapExon_size == 0 && $cumul_overlapIntron_size >0 ){
					
					$partner{'subtype'}		= "intronic";
					$partner{'oversize'}	= $cumul_overlapIntron_size;					
				}				
# 				printf ("Intron :\t%s\t%s\t%d\t%.2f\t%.2f\n", $tr1,$tr2, $cumul_overlapIntron_size, $fractionIntronOver1, $fractionIntronOver2);

								
    	    	# **** CONTAIN ***
				# Compute the overlap size for Contain
				my $cumul_overlapContain_size	=	ExtractFromFeature::intersectExonsOutsideTx( $h1{$tr1}->{'feature'}, $h2{$tr2}->{'feature'}, $s1, $s2, $strandedmode, $verbosity);
				my $fractionContainOver1	=	$cumul_overlapContain_size/$tr1_size;
				my $fractionContainOver2	=	$cumul_overlapContain_size/$tr2_size;				

				if ($cumul_overlapExon_size == 0 && $cumul_overlapIntron_size == 0 && $cumul_overlapContain_size > 0){
					
					$partner{'subtype'}		= "containing";
					$partner{'oversize'}	= $cumul_overlapContain_size;					
				}				

# 				printf ("Contain:\t%s\t%s\t%d\t%.2f\t%.2f\n", $tr1,$tr2, $cumul_overlapContain_size, $fractionContainOver1, $fractionContainOver2);
				
				
			# if intergenic
			} else {

 				$partner{'type'}		= 0 ;					
 				$partner{'dist'}		= Utils::getdistance($b1,$e1,$b2,$e2);					
				$partner{'subtype'}		= intergenic_Subtype($s1, $s2, $b1, $b2);
				$partner{'oversize'}	= 0;
				$partner{'orient'}		= orientation($s1,$s2);
			
			}
			
			push (@{$h_inter{$tr1}->{"partners"}}, \%partner);	
       }
   } 
   return %h_inter;
}

=cut
 Title   : interaction_type
 Usage   : interaction_type($start1, $end1, $start2, $end2)
 Function: define the interation between 2 tx 
 Returns : 1 for genic and 0 for intergenic
 Args    : ranges

=cut
sub interaction_Type{
    my ($b1,$e1,$b2,$e2,$s1,$s2, $stranded) = @_;
	
	if ( Utils::foverlap ($b1,$e1,$b2,$e2, $s1, $s2, $stranded) ){
		return 1;
	} else{
		return 0;
	}
}

=cut
 Title   : orientation
 Usage   : orientation($s1,$2)
 Function: define the orientation between 2 txs
 Returns : sense or antisense
 Args    : strand string

=cut
sub orientation{
    my ($s1,$s2) = @_;
	
	if ($s1 eq "." || $s2 eq "."){ return "NA"}
	return ($s1 eq $s2) ? "sense": "antisense";

}


=cut
 Title   : intergenic_Orientation
 Usage   : intergenic_Orientation($s1, $2, $pos1, $pos2)
 Function: define the orientation between 2 intergenic txs
 Returns : 1 for sense and 0 for antisense
 Args    : strand strings

=cut
sub intergenic_Subtype{
    my ($s1, $s2, $pos1, $pos2) = @_;
	
	die "intergenic_Orientation : feature A position not integer\n" unless ( $pos1 =~ /^-?\d+$/);
	die "intergenic_Orientation : feature B position not integer\n" unless ( $pos2 =~ /^-?\d+$/); 

	if ($s1 eq "." || $s2 eq "."){ return "NA"}

	# if $strand defined
	if ($s1 eq $s2){return "same"}
	else{
		if ($s1 eq "+"){
			if ($pos2 < $pos1){return "div"} else {return "conv"}
		}else{
			if ($pos2 < $pos1){return "conv"} else {return "div"}		
		}
	}
}

1;
