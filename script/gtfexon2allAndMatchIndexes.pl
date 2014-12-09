#!/usr/bin/perl -w


#####################################################################
# TD, Sept 2014
# tderrien@univ-rennes1.fr
# 
# Aims :
#	- Format a gtf to with exon lines to gene and trascript levels
#	- add gene information thanks to Mathieu indexes in /home/genouest/umr6061/recomgen/tderrien/canFam3/annotation/Correspondence_Indexes/
#####################################################################

# Uses
use strict;
use Pod::Usage;
use Getopt::Long;
use Parser;
use Data::Dumper;
use ExtractFromFeature;

# Global Variables
my $infile;
my $mode		= 'all'; # print all attributes for tx and exons
my $help		= 0;
my $verbosity	= 0;
my $man;

# Parsing parameters
my $result = GetOptions(
	"infile|i=s"		=> \$infile,
	"verbosity|v=i"	=> \$verbosity,
	"man|m"			=> \$man,	
	"help|h"		=> \$help);
	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;
pod2usage("I need a gtf file...\n") unless ($infile && -r $infile);
################################################################################
# Indexes
my %RLOC;
my $RLOC_index='/home/genouest/umr6061/recomgen/tderrien/canFam3/annotation/Correspondence_Indexes/RLOCs_index.txt';
open RLOCINDEX ,"$RLOC_index" or die "Error! Cannot open RLOC_index ". $RLOC_index . ": ".$!;
my @lines	=	<RLOCINDEX>;

foreach my $line (@lines){
	chomp ($line);
	my ($rloc, $enscaf, $orth, $biotype) =split (/\t/, $line);
	$RLOC{$rloc}->{'enscaf'}	= $enscaf;
	$RLOC{$rloc}->{'orth'}		= $orth;
}
# ENCAF
my %ENSCAF;
my $ENSCAF_index='/home/genouest/umr6061/recomgen/tderrien/canFam3/annotation/Correspondence_Indexes/ENSCAFGs_index.txt';
open ENSCAFINDEX ,"$ENSCAF_index" or die "Error! Cannot open ENSCAF_index ". $ENSCAF_index . ": ".$!;
@lines	=	<ENSCAFINDEX>;

foreach my $line (@lines){
	chomp ($line);
	my ($enscaf, $ensname, $broadname, $consname, $rloc) =split (/\t/, $line);
	$ENSCAF{$enscaf}->{'ensname'}	= $ensname;
	$ENSCAF{$enscaf}->{'broadname'}	= $broadname;
	$ENSCAF{$enscaf}->{'consname'}	= $consname;
	$ENSCAF{$enscaf}->{'rloc'}		= $rloc;
}



################################################################################
# Parse GTF File at the gene level
print STDERR "Parse Gtf file: '$infile'\n" if ($verbosity > 0);
my $refh              = Parser::parseGTFgene($infile, 'exon',  0, undef, $verbosity);
# print Dumper $refh;

# die;

# Parse gtfHash to be printed
# Gene level
for my $gn (keys %{$refh}){

	# GENE get extra infos
	my $nbtx = scalar (keys %{$refh->{$gn}->{'transcript_id'}});
	
	# biotype list
	my @biotlist;
	foreach my $tr (keys %{$refh->{$gn}->{'transcript_id'}}) {
		push @biotlist , $refh->{$gn}->{'transcript_id'}->{$tr}->{'feature'}[0]->{'transcript_biotype'};
	}
	my $biotlist = join "," , @biotlist;
	
	# Indexes
	my $genename;
	my @genename;
	my $enscaf;	
	my @enscaf;	

	if (exists ($RLOC{$gn}->{'enscaf'})) {
		if ($RLOC{$gn}->{'enscaf'} eq "NA"){
			push @genename, $RLOC{$gn}->{'orth'}; # basically the BROAD mapping of human gene name
			push @enscaf, 'NA';
		} else {
			@enscaf = split ',', $RLOC{$gn}->{'enscaf'}; # if mutliplie enscaf
			foreach my $ensc (@enscaf){
				push @genename, $ENSCAF{$ensc}->{'consname'}; # get all consensus names
			}
		}
	} else {
		die "Error: RLOC $gn does not exist in RLOC_index ". $RLOC_index . "\n";
	}
	$genename	= join ",", @genename;
	$enscaf		= join ",", @enscaf;	
	## Gene level
	print join("\t",$refh->{$gn}->{"chr"}, $refh->{$gn}->{"source"}, "gene", $refh->{$gn}->{"startg"}, $refh->{$gn}->{"endg"}, $refh->{$gn}->{"score"}, $refh->{$gn}->{"strand"}, ".");
 	print "\tgene_id \"".$gn."\"; transcript_nb \"$nbtx\"; transcript_allbiotypes \"$biotlist\"; ensembl_name \"$enscaf\"; gene_name \"$genename\";\n";
	

	## Transcript level	
	foreach my $tr (keys %{$refh->{$gn}->{'transcript_id'}}) {
	
		# shortcut
		my $reftx 		= $refh->{$gn}->{'transcript_id'}->{$tr};
		
		# TX get extra infos
		my $nbex = scalar (@{$reftx->{"feature"}});	
		
	
		print join("\t",$reftx->{"chr"}, $reftx->{"source"}, "transcript", $reftx->{"startt"}, $reftx->{"endt"}, $reftx->{"score"}, $reftx->{"strand"}, ".");
		print "\tgene_id \"".$gn."\"; transcript_id \"".$tr."\"; exon_nb \"$nbex\";";

		# Tx Attrib
		if ($mode eq "all"){
			my %tmph = %{$reftx->{"feature"}[0]}; # take the first exon flags as the flags for the transcripts
			
			# delete unnecesserary keys from hash %tmph
			delete @tmph{qw/feat_level start end strand frame/};
			for (sort keys %tmph){
				print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); 
			}
		}	
		print "\n";
				
		my $countex =0;		
		## Exon level
		foreach my $feat1 (@{$reftx->{"feature"}}) {
		
			$countex++;
			
			print join("\t",$reftx->{"chr"}, $reftx->{"source"}, "exon", $feat1->{"start"}, $feat1->{"end"}, $reftx->{"score"}, $reftx->{"strand"}, $feat1->{"frame"});
			print "\tgene_id \"".$gn."\"; transcript_id \"".$tr."\"; exon_number \"$countex\";";
 			
			# Ex Attrib
			if ($mode eq "all"){
				my %tmph = %{$feat1};
				# delete unnecesserary keys from hash %tmph
				delete @tmph{qw/feat_level start end strand frame/};
				for (sort keys %tmph){
					print " $_ \"$tmph{$_}\";" if (defined $tmph{$_});
				}
			}	
			print "\n";
		}		
	}
}

__END__

=head1 NAME

gtfexon2allAndMatchIndexes.pl - Format a .gtf file with exon levels  in 3 levels exon, transcript and gene and add matching gene with Indexes

=head1 SYNOPSIS

perl gtfexon2allAndMatchIndexes.pl -infile <gtf file> [Options] 

Format a .gtf file with exon levels  in 3 levels exon, transcript and gene

Options:

	-help		: Help message

	-man		: man help

	-verbosity	: level of verbosity 

	
=head1 OPTIONS

=over 8

=item B<-verbosity>
: Level of verbosity to follow process

=item B<-man>
: Print man help and exit

=item B<-help>
: Print help message and exit

=back

=head1 DESCRIPTION

=over 8

=item * Parse a gtf file with *only* exon levels

=item * Format the file with exon, transcript and gene levels in gtf

=item * Parse Indexes and add correspondent gene name and enscaf

=back

=cut
