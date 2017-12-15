#!/usr/bin/perl -w
#
#
#######################################################
# Author  :  Thomas DERRIEN 
# email   :  tderrien@univ-rennes1.fr
# Purpose :  Compute specificity score as described in 
#				* in  Yanai I. et al., Bioinformatics , 2005, vol. 21 (pg. 650-659)
#				* see also : Kryuchkova-Mostacci, N. & Robinson-Rechavi, M. A benchmark of gene expression tissue-specificity metrics. Brief. Bioinform. 18, 205–214 (2017).

########
#
#
# Usage : 
# 	 specificity_score.pl matrix.rpkm | column -t 
#
# Example matrix.rpkm
# ID           HAIRFOLLICULE  ADRENALGLAND  CORTEX  GUTCOLON  JEJUNUM  KERATINOCYTE  MAMMARYGLAND_GSMD  OLFBULB_GSMD  PANCREAS
# XLOC_000019  47             552           91      15        93       0             29                 14            6
# XLOC_000030  46             20            10      22        3        12            20                 0             1
# XLOC_000031  0              239           0       0         0        0             0                  0             0
# XLOC_000032  1              0             57      0         0        0             0                  190           0
# XLOC_000041  218            301           415     350       171      67            366                246           71
# XLOC_000052  0              0             60      7         2        0             0                  169           0
########################################################


use strict;
use warnings;
use Data::Dumper;
use List::Util qw(max sum first);

my $infile 	= shift or die "Usage: $0 MATRIX_FILE\n";
my $minexpr =   1; # minimum expression level 

open (my $fh , "<" , $infile) or die "Cannot open < $infile: $!";

my $header=<$fh>;   #First line is read here 

# get the first line in an array
my @header = split /\s+/,$header;

while (<$fh>){

	my @ar		=	split /\s+/;
	my $id 		= 	$ar[0];
	my @vec		=	@ar[1 .. $#ar]; # get a slice until the last element of the array
	
    my $maxexp	=	max(@vec);
    
    next if ($maxexp < $minexpr);

	# Compute the scpeficity score
	my ($sp_score, $tissue)    =   getTSscore (\@vec, \@header, $minexpr);
	
	print "$id\t$sp_score\t$tissue\t", join("\t",@vec),"\n";
	
}

close ($fh);


# Subfunction to compute as in Yanai et al. : sum(1-(d/max(d)))/(length(d)-1)
sub getTSscore {

    my ($vec, $header)  =   @_;
    $minexpr    ||=  0;
    
    # output variables
    my $sp_score    = 0;
    my $besttissue  = "NA";    
    
    # deref array ref
    my @vec     =   @{$vec};
    my @header  =   @{$header};
    
    # get max per tissue
	my $maxexp	=	max(@vec);
	
	# test if 0 vector
	return $sp_score, $besttissue if ($maxexp == 0);

	my $nbexp	= 	@vec -1; 		# size - 1
	my $sum		=	0;
	foreach my $exp (@vec){
		$sum += 1 - ($exp/$maxexp);
	}
	# divide by the number of exp - 1
	$sp_score	= $sum /$nbexp;    

	# get indice of the max tissue
	my $maxind	    =	first { $vec[$_] eq $maxexp } 0..$#vec;
	$besttissue  =   $header-> [$maxind+1];
	
    return $sp_score, $besttissue;
}


