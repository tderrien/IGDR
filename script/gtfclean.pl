#!/usr/bin/perl -w


# Aim :clean a gtf file
###############################################################

=head1 NAME

gtfclean.pl : clean a gtf file for duplicate transcripts (same exons), matched transcripts (same intron chains), remove small introns...

=head1 SYNOPSIS

gtfclean.pl  -i infile.gtf [options]

=head1 OPTIONS

=over 8

=item B<-i|infile>

Input .gtf file. [mandatory]

=item B<-ke|keepduptx>

keep duplicate transcripts having same exons... [default: FALSE (only one representative isoform is kept)]

=item B<-ki|keepdupint>

keep duplicate transcripts having same introns chain... [default: TRUE (use --no-keepdupint to keep only on representative isoform)]

=item B<-m|merge_intron>

Introns smaller than this value will be merged...  [ default: 5 nt].

=item B<-p|proc>

number of process [default 1].

=item B<-l|logfile>

Write the result of the cleaning in prefix.log file [default: STDERR]

=item B<-v|verbosity>

Level of verbosity [default 0].

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will parse a .gtf file and return a cleaned one

=cut


# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;
use File::Basename;


# Own lib : /home/genouest/umr6061/recomgen/tderrien/bin/ThomasPerl/lib/
use Parser;
use ExtractFromHash;
use Intersect;
use Utils;

use Parallel::ForkManager;


my $progname=basename($0);

# Variables
my $infile		='';
my $keepduptx 	= 0 ;
my $keepdupint 	= 1;
my $proc 		= 4;
my $merge_intron= 5;
my $logfile		= 0;
my $man 		= 0;
my $help 		= 0;
my $verbosity	= 0;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
	'i|infile=s'	 	=> \$infile,
   'ke|keepduptx!' 		=> \$keepduptx,
   'ki|keepdupint!' 	=> \$keepdupint,
   'm|merge_intron=i' 	=> \$merge_intron,
   'p|proc=i' 			=> \$proc,
   'l|logfile=s' 		=> \$logfile,
   'v|verbosity=i'		=> \$verbosity,
   'help|?' 			=> \$help,
   'man' 				=> \$man
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

# print STDERR "Options:\n- infile: $infile\n- keepduptx: $keepduptx\n- keepdupint: $keepdupint\n- merge_intron: $merge_intron\n";

# Test parameters
pod2usage("Error: Cannot read your input .gtf file '$infile'...\nFor help, see:\n$progname --help\n") unless( -r $infile);


# Define log file name
my $logfilename;
my $fh;
if ($logfile){
	if ($logfile eq ""){
		my $basename = basename($infile);
		($logfilename = $basename) =~ s/\.[^.]+$/.log/;
	} else {
		$logfilename = $logfile;
	}
	# open fh
	open($fh, '>', $logfilename) or die "Could not open file '$logfilename' $!";
}


my $splitbychr	=	1;
# Parsing
# my %hA_chr  = Parser::parseGTFSplitchr($infile, 'exon', $verbosity);
my $hA_chr		= Parser::parseGTF($infile, 'exon',  1,1 , $verbosity);

my $hB_chr 	= $hA_chr;

# Launch parralel
my $pm = Parallel::ForkManager->new($proc);
my %good_tx = ();  # for collecting responses
$pm -> run_on_finish (
  sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result_ref) = @_;
		my $chr = $ident;

		warn("Child $chr killed by signal $exit_signal"), return if $exit_signal;
		warn("Child $chr exited with error $exit_code"),  return if $exit_code;
		warn("Child $chr encountered an unknown error"),  return if !$result_ref;

	   $good_tx{$chr} = $result_ref;
  }
);

				
# Process 2 hashes
foreach my $chrA (keys %{$hA_chr} ) {

	my %h;
	
	foreach my $chrB (keys %{$hB_chr} ) {
	
		next if ($chrA ne $chrB);
                
		# start fork
	    my $pid = $pm->start($chrB) and next;
     	
     	# get ref on hash per chromosome
     	my $refh1 = $hA_chr->{$chrA};
        my $refh2 = $hB_chr->{$chrB};        
        
        # Launch interactionCollection
		print STDERR "$chrB\n" if ($verbosity > 5);
		# Launch Intersect
        %h = Intersect::Intersect2HsplitFork_clean($refh1, $refh2, $keepduptx, $keepdupint, $fh, $verbosity);

    }
    # end fork
    $pm->finish(0, \%h);

}


# wait all sub process
$pm->wait_all_children;


# Parse good transcripts by chr
foreach my $chr (keys %good_tx ) {

	my %chrh = ExtractFromHash::MergeSmallIntrons($good_tx{$chr}, $merge_intron);
	ExtractFromHash::printGTF(\%chrh, 'all', $verbosity);
}

if (defined $fh){
	close $fh;
}