#!/usr/bin/perl -w


#####################################################################
# TD, Sept 2020
# tderrien@univ-rennes1.fr
# 
# Aims :
#	- Format a gtf with missing tx/ex strand (e.g. - or *), attribute gene strand to tx/ex strand or "+" default
#####################################################################

# Uses
use strict;
use Pod::Usage;
use Getopt::Long;
use Parser;
use Data::Dumper;


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
# Parse GTF File at the gene level
print STDERR "Parse Gtf file: '$infile'\n" if ($verbosity > 0);
my $refh              = Parser::parseGTFgene($infile, 'exon',  0, undef, undef, $verbosity);

# Parse gtfHash to be printed
# Gene level
for my $gn (keys %{$refh}){

	## Gene level
#	print join("\t",$refh->{$gn}->{"chr"}, $refh->{$gn}->{"source"}, "gene", $refh->{$gn}->{"startg"}, $refh->{$gn}->{"endg"}, $refh->{$gn}->{"score"}, $refh->{$gn}->{"strand"}, ".");
# 	print "\tgene_id \"".$gn."\";\n";
		

	## Transcript level	
	foreach my $tr (keys %{$refh->{$gn}->{'transcript_id'}}) {
	
		# shortcut
		my $reftx 		= $refh->{$gn}->{'transcript_id'}->{$tr};
		my $flag		= "restrand_ok";
	
		if ($reftx->{"strand"} eq "." || $reftx->{"strand"} eq "*") {
			if ( $refh->{$gn}->{"strand"} ne "."){ 
				$reftx->{"strand"} = $refh->{$gn}->{"strand"}; # tx restrand is = to gene restrand
				$flag              = "restrand_gn"
			} else {
				$reftx->{"strand"} = "+"; # default restrand
				$flag 		   = "restrand_+"
			}
		}
	
		print join("\t",$reftx->{"chr"}, $reftx->{"source"}, "transcript", $reftx->{"startt"}, $reftx->{"endt"}, $reftx->{"score"}, $reftx->{"strand"}, ".");
		print "\tgene_id \"".$gn."\"; transcript_id \"".$tr."\";";

		# Tx Attrib
		if ($mode eq "all"){
			my %tmph = %{$reftx->{"feature"}[0]}; # take the first exon flags as the flags for the transcripts
			# delete unnecesserary keys from hash %tmph
			delete @tmph{qw/feat_level start end strand frame/};
			for (sort keys %tmph){
				print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); 
			}
		}	
		print " flag_st \"".$flag."\";";
		print "\n";
						
		## Exon level
		foreach my $feat1 (@{$reftx->{"feature"}}) {
		
			print join("\t",$reftx->{"chr"}, $reftx->{"source"}, "exon", $feat1->{"start"}, $feat1->{"end"}, $reftx->{"score"}, $reftx->{"strand"}, $feat1->{"frame"});
			print "\tgene_id \"".$gn."\"; transcript_id \"".$tr."\";";
 			
			# Ex Attrib
			if ($mode eq "all"){
				my %tmph = %{$feat1};
				# delete unnecesserary keys from hash %tmph
				delete @tmph{qw/feat_level start end strand frame/};
				for (sort keys %tmph){
					print " $_ \"$tmph{$_}\";" if (defined $tmph{$_});
				}
			}	
                        print " flag_st \"".$flag."\";";
			print "\n";
		}		
	}
}

__END__

=head1 NAME

gtf2restrand.pl - Format a .gtf file with strand = "*" or ".", attribute gene strand to tx/ex or "+" (default)

=head1 SYNOPSIS

perl gtf2restrand.pl  -infile <gtf file> [Options] 

Format a .gtf file with strand = "*" or ".", attribute gene strand to tx/ex or "+" (default)

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

=item * Format a .gtf file with strand = "*" or ".", attribute gene strand to tx/ex or "+" (default)

=back

=cut
