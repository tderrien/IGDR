#!/usr/bin/perl -w



########################################################
# TD, July 2013
# tderrien@univ-rennes1.fr
# 
# Aims :
#	- parse a gtf file A
#	- parse a double gtf file B from intersectRetouch.pl
# 	- add matching tx attributes (matchEns75Tx100 matchEns75Gn100 from) from B
#########################################################

# Uses
use strict;
use Pod::Usage;
use Getopt::Long;
use File::Basename;


# Own package
use StringUtils;
use Utils;
use Parser;
use ExtractFromFeature;
use Data::Dumper;

# Variables
my $infileA		='';
my $infileB		='';

my $help 		= 0;
my $man 		= 0;
my $verbosity	= 0;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
	'a|infileA=s'	 	=> \$infileA,
	'b|infileB=s'	 	=> \$infileB,	
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'man' 				=> \$man
) or pod2usage(2);	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-exitval=>0, -verbose=>2) if $man;

my $progname=basename($0);


# Test parameters
pod2usage("Error: Cannot read your input A .gtf file '$infileA'...\nFor help, see:\n$progname --help\n") unless( -r $infileA);
pod2usage("Error: Cannot read your input B .gtf file '$infileB'...\nFor help, see:\n$progname --help\n") unless( -r $infileB);



my $splitbychr	=	0;
my $parseextraf	=	undef;

# Parsing
my $hA_chr		= Parser::parseGTF($infileA, 'exon',  $splitbychr, $parseextraf , $verbosity);
my $hB_chr		= Parser::parsedoubleGTF($infileB, 'exon',  0, undef, $verbosity);


##################################################################################
# Get Sequences
my $i=0;
my $sizeh	= keys(%{$hA_chr}); 


for my $tr (keys %{$hA_chr}) {
	
	# if matching tx from B
	my $refh;
	if (exists($hB_chr->{$tr})){
		print STDERR "Exist $tr\n" if ($verbosity >5);
		$refh = ExtractFromFeature::addKeyValAttrib($hA_chr->{$tr},  $hB_chr->{$tr}, 'matchEns75Tx100', 'matchEns75Gn100', $verbosity)
	} else{
		print STDERR "Does not Exist $tr\n" if ($verbosity >5);
		$refh = ExtractFromFeature::addKeyValAttrib($hA_chr->{$tr},  0, 'matchEns75Tx100', 'matchEns75Gn100', $verbosity)
	}

	# assign to new hash to be printed
	$hA_chr->{$tr}= $refh;
	
	if ($verbosity > 0){
		Utils::showProgress( $sizeh , $i++, "AddKeyValAttrib :");
	}				
}

# print Dumper $hA_chr;
 ExtractFromHash::printGTF($hA_chr, 'all', $verbosity);



=head1 NAME

addAttribfromIntersect.pl :  add matching tx attributes (matchEns75Tx100 matchEns75Gn100from) from B to A 

=head1 SYNOPSIS

addAttribfromIntersect.pl  -a infileA.gtf -b infileB.gtf [options]

=head1 OPTIONS

=over 8

=item B<-a|afile>

Input A .gtf file. [mandatory]

=item B<-b|bfile>

Input B double .gtf file from intersectRetouch. [mandatory]

=cut