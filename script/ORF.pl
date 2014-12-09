#!/usr/bin/perl -w


use strict;
use warnings;

# general lib
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Bio::DB::Fasta;

# Lib perso
use Parser;
use StringUtils;
use Utils;
use ExtractFromFeature;
use ExtractFromHash;
use Orf;

# Variables
my $infile		= "";
my $genome		= "";
my $tx_biotype	= "";
my $man 		= 0;
my $help 		= 0;
my $verbosity	= 0;

my $progname=basename($0);


# Prog var
my $splitbychr = 0;


####################
# GetOptions 
GetOptions(
	'i|infile=s'	 	=> \$infile,
	'g|genome=s' 		=> \$genome,
	'b|biotype:s' 		=> \$tx_biotype,
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'man' 				=> \$man
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test parameters
pod2usage("Error: Cannot read your input .gtf file '$infile'...\nFor help, see:\n$progname --help\n") unless( -r $infile);
pod2usage("Error: Cannot read your genome file '$genome'...\nFor help, see:\n$progname --help\n") if (! -r $genome && !-d $genome);


#######################
# Parse GTF file
my $h		= Parser::parseGTF($infile, 'exon', $splitbychr, $tx_biotype , $verbosity);

my $sizeh = keys(%{$h});
my $i=0;

for my $tr (keys(%{$h})){

	my $chr = $h->{$tr}->{'chr'};
	my $strand = $h->{$tr}->{'strand'};

	# get cDNA sequence for transcript tr
	my $seq = ExtractFromFeature::feature2seq($h->{$tr}->{'feature'}, $genome, $chr , $strand, 0, $verbosity);
	
	# select best ORF
	my $orf_selected =  Orf::chooseORF($seq);
	
	# get values of the transcript with CDS/UTR
	my $refh = ExtractFromFeature::Tx2CDS($h->{$tr}, $orf_selected);

	# assign to new hash to be printed
	$h->{$tr}= $refh;
	
	if ($verbosity > 0){
		Utils::showProgress( $sizeh , $i++, "Compute ORF and Project:");
	}				
}

ExtractFromHash::printGTF($h, 'all', $verbosity);


__END__

=head1 NAME

ORF.pl : compute the longest ORF of transcripts in .gtf file and using a genome fasta file

=head1 SYNOPSIS

ORF.pl  -i infile.gtf -g genome.fa [options]

=head1 OPTIONS

=over 8

=item B<-i|infile>

Input .gtf file. [mandatory]

=item B<-g|genome>

Input genome file. [mandatory]

full path to a multi-fasta file with the genomic sequences
      for all input features, OR a directory with single-fasta files
      (one per genomic sequence, with file names matching sequence names)


=item B<-v|verbosity>

Level of verbosity [default 0].

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> parse a .gtf file and return the longest ORF for each transcript

=cut
