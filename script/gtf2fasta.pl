#!/usr/bin/perl -w



########################################################
# TD, May 2013
# tderrien@univ-rennes1.fr
# 
# Aims :
#	- Parse a gtf files
#	- Get Fasta sequences of concatened exons
# Note
#	- Updated version of -sav with personal perl packages in : /home/genouest/umr6061/recomgen/tderrien/bin/ThomasPerl/lib/
#########################################################

# Uses
use strict;
use Pod::Usage;
use Getopt::Long;

#use lib '/home/genouest/umr6061/recomgen/tderrien/bin/perl/lib/';
# use lib '/home/genouest/umr6061/recomgen/tderrien/src/bioperl-1.2.3';

# Own package
use StringUtils;
use Utils;
use Parser;
use ExtractFromFeature;
use Bio::DB::Fasta;
use Data::Dumper;
use Bio::SeqIO;

# Global Variables
my $infile;
my $genome    ="/home/genouest/umr6061/recomgen/tderrien/canFam3/sequence/softmasked/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa";
# my $genome    ="/omaha-beach/tderrien/DATA/hg19/"; # for human
my $slop                = 0;        # nb of extranucleotide around transcript
my $outfile             = undef;    # name of the output file (defaut stdout)
my $cds					= 0;		# only extract CDS sequence (from start to stop) [default FALSE]
my $fullcds				= 0;		# only extract full length CDS sequence i.e only those that have a START and END codon [default FALSE]

my $help                = 0;
my $verbose             = 0;
my $man;

my %filtertag;


# Parsing parameters
GetOptions(
	"infile|i=s"		=> \$infile,
	"genome|g=s"		=> \$genome,
	"outfile|o=s"		=> \$outfile,
	"slop|s=i"			=> \$slop,
	"f|filter=s"        => \%filtertag,	
	"cds|c!"			=> \$cds,	
	"fullcds!"			=> \$fullcds,	
	"verbose|v=i"		=> \$verbose,	
	"man|m"				=> \$man,	
	"help|?"			=> \$help) or pod2usage(2);	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-exitval=>0, -verbose=>2) if $man;

# Mandatory options 
pod2usage("I need a gtf file") unless ($infile && -r $infile);
pod2usage("The slop value should be a positive number and lower than 1Mb") unless ($slop >=0 && $slop <1000000);
pod2usage("I need a valid genome file or directory: '$genome'\n") if  (! -r $genome && ! -d $genome);
################################################################################


# Parse file at the exon level
my $h;		
if ($cds || $fullcds){
	warn "CDS option 'cds=$cds , fullcds=$fullcds' activated...\n **only CDS and stop_codon lines will be parsed**\n";
	sleep(3);
	$h	=	Parser::parseGTF($infile, 'CDS,stop_codon', 0, \%filtertag , $verbose); # WARNING, we still include stop codon since it is not part of the CDS by defintion
}else{
	$h	= Parser::parseGTF($infile, 'exon', 0, \%filtertag , $verbose);
}



##################################################################################
# Get Sequences
my $i=0;
my $h_transcript_size	= keys(%{$h}); 
my $seqdata;

if ($verbose >0){ print STDERR "Get Sequences\n";}

# create genome DB
my $db       = Bio::DB::Fasta->new($genome);

my $seqOUT;
# Output
if (defined $outfile){
	$seqOUT	=	Bio::SeqIO ->new(-format => 'fasta', -file => ">$outfile", -alphabet =>'dna');
} else {
	$seqOUT	=	Bio::SeqIO ->new(-format => 'fasta', -fh     => \*STDOUT, -alphabet =>'dna');
}

for my $tr (keys %{$h}){
	
	#Initalize sequence:
	my $seqstring	=	"";
	my $id_sequence	=	"";
	my $cpt			= 	0;
	my $chr 		=	$h->{$tr}->{"chr"};
	
	if (!defined($db->seq($chr, 0 => 10))){
		warn "WARNING: $tr is on chr '$chr' which is not a chromosome of the genome fasta file : $genome... Skip this tx!\n";
		next;
	}

	
	foreach my $exon (@{$h->{$tr}->{"feature"}}) {
		
		$cpt++;
		my $s	= $exon->{"start"};
		my $e	= $exon->{"end"};
		
		if ( $cpt == 1){ 												# if first, we add $slop bp to START
			$seqstring   	.=  $db->seq($chr, $s - $slop => $e);
		} elsif ($cpt == scalar(@{$h->{$tr}->{"feature"}})) {
			$seqstring   	.=  $db->seq($chr, $s => $e+$slop); 		# if last exon, we add $slop bp to END
		} else{
			$seqstring   	.=  $db->seq($chr, $s => $e);				# else no slop
		}
    }
    
    
	#RevComp if strand -
	if ( $h->{$tr}->{"strand"} eq '-'|| $h->{$tr}->{"strand"} eq '-1') {
        $seqstring 	= StringUtils::getRevComp($seqstring);
	}


	# Summarize data e.g >TCONS_00005869 XLOC_001028_-_1:2753268-2784339_Cufflinks
	# and fasta sequence
	# header
	my $tx_biot	=	ExtractFromFeature::getKeyFromFeature($h, $tr, 'transcript_biotype', $verbose);
	my $id		= $tr."_".$h->{$tr}->{"gene_id"}."_".$h->{$tr}->{"strand"}."_".$h->{$tr}->{"chr"}.":".$h->{$tr}->{"startt"}."-".$h->{$tr}->{"endt"}."_".$tx_biot;
	my $new_seq = Bio::Seq->new(-id => $id, -seq => $seqstring);
	
	# If full CDS has to be checked

	if ($fullcds){
		
		if (checkDnaStartStop($new_seq->seq())){
			$seqOUT->write_seq($new_seq);
		} else{
		
			warn "WARNINGS: $tr: Option --fullcds is activated but start and/or stop codons is not found\n";
		}
	} else {
		$seqOUT->write_seq($new_seq);		
	}

	if ($verbose > 0){
        Utils::showProgress($h_transcript_size, $i++, "Print ".$tr.": ");
    }

}


# Check if full CDS
sub checkDnaStartStop {
my ($seq) = @_;
	
	if (defined ($seq) && $seq=~m/^ATG/i && $seq =~m/TAA$|TAG$|TGA$/i){
		return 1;
	} else {
		return 0;
	}
}

__END__

=head1 NAME

gtf2fasta.pl - Extract fasta sequence from a gtf file


=head1 SYNOPSIS

perl gtf2fasta.pl -infile <gtf file> [Options]



Options:

	--genome|-g	: path directory to genome sequence file/dir  [default: ~/data/canFam3/sequence/softmasked/Canis_familiaris.CanFam3.1.72.dna_sm.toplevel.fa]

	--filter|-f	: filtering attributes such as key=val (-f chr=X,34 -f transcript_biotype=lincRNA,antisense) [default: undef]

	--slop|-s	: number of nucleotide to be added in 5' and 3' of the transcript sequence [default: 0]

	--cds|-c	: extract CDS sequence wrt to CDS features in the .GTF files (may be incomplete) [default: FALSE]

	--fullcds	: extract *ONLY* full length CDS sequence wrt to CDS, start_codon, stop_codon features in the .GTF files [default: FALSE]
	
	--outfile|-o	: output file name [default : STDOUT]
	
	--help|-h	: Help message [default : 0]

	--man|-m	: man help [default : 0]

	--verbose|v 	: level of verbosity [default : 0]

	
=head1 OPTIONS

=over 8

=item B<--genome|-g>
: Path to the directory that contains fasta sequences of species chromosomes

=item B<--filter|-f>
: filtering attributes such as key/val ( -f transcript_biotype=protein_coding,lincRNA -f transcript_id=ENSTA,ENSTB)

=item B<--outfile|-o>
: output file name

=item B<--cds|-c>
: extract CDS sequence wrt to CDS features in the .GTF files [default: FALSE]

=item B<--fullcds>
: extract *ONLY* full length CDS sequence wrt to CDS, start_codon, stop_codon features in the .GTF files [default: FALSE]

=item B<--slop|-s>
: number of nucleotide to be added in 5' and 3' of the transcript sequence

=item B<--verbose|-v>
: Level of verbosity to follow process

=item B<--man|-m>
: Print man help and exit

=item B<--help|-h>
: Print help message and exit

=back

=head1 DESCRIPTION

=over 8

=item * Parse a gtf files and a genome file/directory using Bio::DB::Fasta

=item * Get Transcripts Fasta sequences of concatened exons

=item * Filter for some particular tags (-f transcript_biotype)

=cut
