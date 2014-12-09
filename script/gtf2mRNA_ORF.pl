#!/usr/bin/perl -w

#
# Step 2014
# tderrien@univ-rennes1.fr

# - Take . gtf file 
# - Compute ORF hash object if
#		- CDS is  given
#		- by longestORF module otherwise
# - write 2 files (cDNA) and ORF if max num of ORF (numtx) is reached 

# ToDo
# Add a min size ORF if longestORF module
##########################################################################################


=head1 NAME

gtf2mRNA_ORF.pl - Extract mRNA and ORF (complete) fasta sequences from a .gtf file

=head1 SYNOPSIS

gtf2mRNA_ORF.pl  -i infile.gtf -g genome.fa  [options]

  Options:
	--infile|-i	: input mRNA .gtf file [mandatory]
	--genome|-g	: input genome file (either a file or directory containing chr files) [mandatory]
	--cdna|-c	: name of the output cDNA/mRNA file [default 'cdnafile.fa' ]
	--orf|-o	: name of the output ORF file [default 'orf.fa']
	--numtx|-n	: number of transcripts to extract [default 2000 ]
	--verbose|-v	: level of verbosity [default 0 ]
	--help		: brief help message
	--man|-m	: full documentation

=head1 OPTIONS

=over 8

=item B<-cdna|c>

name of the cDNA/mRNA output fasta file

=item B<-orf|o>

name of the ORF output fasta file

=item B<-ntx|n>

Integer corresponding to the  number of transcript sequences to be extracted by the program for the training

=item B<-verbose|v>

The level of verbosity of the program 

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input .gtf file and extract mRNA and ORF sequences

=cut




############################ library #####################################################
use warnings;
use strict;

use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Bio::SeqIO;
use Bio::DB::Fasta;

use Parser;
use ExtractFromFeature;
use Orf;
use Utils;


# Variables
my $infile		= "";
my $genome		= "";
my $cdnafile	= "cdnafile.fa";
my $orffile		= "orffile.fa";
my $numtx		=	undef; # min nb of complete ORF to reach
my $man 		= 0;
my $help 		= 0;
my $verbosity	= 0;

GetOptions (
	"infile|i=s" 	=>	\$infile,
	"genome|g=s"	=>	\$genome,
	"cdna|c=s"		=>	\$cdnafile,
	"orf|o=s"		=>	\$orffile,		
	"ntx|n=i"		=>	\$numtx,
	"verbose|v=i"	=>	\$verbosity,		
	"help|?"		=>	\$help,
	"man|m"			=>	\$man,
	
) or pod2usage(2);	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-exitval=>0, -verbose=>2) if $man;

my $progname=basename($0);


# Test parameters
pod2usage("Error: Cannot read your input .gtf file '$infile' (-i option)...\nFor help, see:\n$progname --help\n") unless( -r $infile);
pod2usage("Error: Cannot read your genome file '$genome' (-g option)...\nFor help, see:\n$progname --help\n") if (! -r $genome && !-d $genome);


#######################
# Parse GTF file
#   * with CDS
#   * with stop codon for checking ORF integrity
my $h		= Parser::parseGTF($infile, 'exon,CDS,stop_codon', 0, undef , $verbosity);
my $sizeh = keys(%{$h});
die "Your input gtf '$infile' contains only *$sizeh* transcripts.\nNot enough to training the program (default option --ntx|-n)\n" if (defined $numtx && $sizeh < $numtx);

print STDERR "Your input gtf '$infile' contains *$sizeh* transcripts\n" if ($verbosity > 0 );


##############################################
my $orfob;
my $allow_no_start 	=	0;	# do not allow for CDS start not found
my $allow_no_stop 	=	0;	# do not allow for CDS stop not found
my %h_orf;					# for storing and printing ORF sequence
my $countorf		=	0;	# counter on good ORF (start and end found)
my $filterforCDS	=	0;	# get only line with CDS level


# Output files
my $cdnaseqOUT	=	Bio::SeqIO ->new(-format => 'fasta', -file => ">$cdnafile", -alphabet =>'dna');

my $i = 0;
for my $tr (keys(%{$h})){

	print STDERR "Tx:", $i++,"/$sizeh\r" if ($verbosity>5);
	# shortcut for feature2seq sub
	my $chr 	= $h->{$tr}->{'chr'};
	my $strand 	= $h->{$tr}->{'strand'};

	# Check Biotype
	my $biotype = $h->{$tr}->{'feature'}[0]->{'transcript_biotype'} if (defined $h->{$tr}->{'feature'}[0]->{'transcript_biotype'});

 
	# get cDNA sequence for transcript tr
	$filterforCDS		=	0; # do we filter seq for CDS
	my $cdnaseq 		=	ExtractFromFeature::feature2seq($h->{$tr}->{'feature'}, $genome, $chr , $strand, $filterforCDS, $verbosity);	
            
	#######################################
	# ORF
	my $containCDS =   ExtractFromFeature::checkCDS($h->{$tr}->{'feature'});
	if (! $containCDS ){
		# we create an ORF hash 	
		$orfob	=	Orf::longestORF2($cdnaseq,".", $allow_no_start, $allow_no_stop,0,1);
# 		print Dumper $orfob;
# 		print  "$tr ",length($orfob->{'cds_seq'}) -3,"\n";
	} else {
		$filterforCDS	= 1; # we activate filter to get only CDS and stop codon DNA sequence
		my $orfseq		= ExtractFromFeature::feature2seq($h->{$tr}->{'feature'}, $genome, $chr , $strand, $filterforCDS, $verbosity);
		# we create an ORF hash 
		$orfob			= Orf::orfSeq2orfOb($orfseq, $strand, $verbosity);
		
	}
	    
	# Add ORF to a hash %h_orf only if the ORF is complete
	if ($orfob->{'check_start'} && $orfob->{'check_stop'}){
		$h_orf{$tr}	=	$orfob->{'cds_seq'};
 		print STDERR "Extracting ", $countorf++,"/$numtx...\r" if (defined $numtx);
 		    
	} else {
		warn "Tx: $tr ('$biotype') with CDS features: $containCDS is not complete...skipping for training\n" if ($verbosity > 5);
		next; # next if ORF is not OK
	}
	
	#######################################
	# ADD cDNA only if ORF is OK
	# Write cDNA seq
	my $new_seq = Bio::Seq->new(-id => $tr, -seq => $cdnaseq, -alphabet => 'dna');
    $cdnaseqOUT->write_seq($new_seq); 	
	
	if (defined $numtx && $countorf == $numtx){
		print STDERR "Max ORF sequences '$numtx' reached..ending!\n";
		last;
	}
	
	
}

# Final Check if the number of complete ORF is ok
my $sizehorf = keys(%h_orf);
die "The number of complete ORF found with computeORF mode is *$sizehorf* transcripts... That's not enough to training the program\n" if (defined $numtx && $sizeh < $numtx);

print STDERR "Writing orf file '$orffile'\n" if ($verbosity > 5);

my $orfseqOUT	=	Bio::SeqIO ->new(-format => 'fasta', -file => ">$orffile", -alphabet =>'dna');

foreach my $orfid (keys %h_orf){
    my $new_seq = Bio::Seq->new(-id => $orfid, -seq => $h_orf{$orfid});
    $orfseqOUT->write_seq($new_seq);	
}



