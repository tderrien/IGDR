#!/usr/bin/perl -w


# Aim : compare 2 gtf files and return the unique portion
###############################################################

# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;



# Own lib : /home/genouest/umr6061/recomgen/tderrien/bin/ThomasPerl/lib/
use Parser;
use ExtractFromHash;
use Intersect;
use Utils;

use Parallel::ForkManager;


my $progname=basename($0);

# Variables
my $infileA		='';
my $infileB		='';
my $uniquexa 	= 0 ;
my $uniquexb 	= 0 ;
my $common		= 1;
my $fraction 	= 1 ;
my $proc 		= 4;
my $strandedmode= 0; # default unstranded
my $man 		= 0;
my $help 		= 0;
my $verbosity	= 0;

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
	'a|infileA=s'	 	=> \$infileA,
	'b|infileB=s'	 	=> \$infileB,	
	'ea|uniquexa!' 		=> \$uniquexa,
	'eb|uniquexb!' 		=> \$uniquexb,
	'c|common!' 		=> \$common,
	'f|fraction=f' 		=> \$fraction,	
	's|strandmode=i' 	=> \$strandedmode,
	'p|proc=i' 			=> \$proc,
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'man' 				=> \$man
) or pod2usage(2);	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-exitval=>0, -verbose=>2) if $man;

print STDERR "Options:\n- infileA: $infileA\n- infileB: $infileB\n- uniquexa: $uniquexa\n- uniquexb: $uniquexb\n- common: $common\n" if ($verbosity > 5);

# Test parameters
pod2usage("Error: Cannot read your input A .gtf file '$infileA'...\nFor help, see:\n$progname --help\n") unless( -r $infileA);
pod2usage("Error: Cannot read your input B .gtf file '$infileB'...\nFor help, see:\n$progname --help\n") unless( -r $infileB);
pod2usage ("- Error: \$fraction option '$fraction' should be a float between 0 and 1 [0-1] (e.g 0.5 if 50% overlap is required)\n") unless ($fraction >= 0 and $fraction <= 1);

pod2usage("Error: Cannot report common tx while unique A or B activated. Use '--no-common' option to exclude matching tx.\n\nFor help, see:\n$progname --help\n") if( ($common && $uniquexa) || ($common && $uniquexb));
pod2usage ("- Error: \$strandedmode option '$strandedmode' should be 0 (unstranded mode comparison), 1 (same strand) or -1 (for opposite strand e.g antisense)\n") if ($strandedmode !=0 && $strandedmode !=1 && $strandedmode != -1 );


my $splitbychr	=	1;
my $parseextraf	=	undef;

# Parsing
my $hA_chr		= Parser::parseGTF($infileA, 'exon',  $splitbychr, $parseextraf , $verbosity);
my $hB_chr		= Parser::parseGTF($infileB, 'exon',  $splitbychr, $parseextraf , $verbosity);


# Launch parralel
my $pm = Parallel::ForkManager->new($proc);
my %retrieved_responses = ();  # for collecting responses
$pm -> run_on_finish (
  sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result_ref) = @_;
		my $chr = $ident;

		warn("Child $chr killed by signal $exit_signal"), return if $exit_signal;
		warn("Child $chr exited with error $exit_code"),  return if $exit_code;
		warn("Child $chr encountered an unknown error"),  return if !$result_ref;

	   $retrieved_responses{$chr} = $result_ref;
  }
);

# For option uniqueA or B, we need to check whether some txs are specific of ONE chromosome
# that will be miss by the chr to chr comparison
my $refarray_chrA;
my $refarray_chrB;	
if (!$common){
	# get chromosomes of A that are not in B in an arrayref and reciproc
	if ($uniquexa){
		$refarray_chrA	=	compare_chr_distrib($hA_chr, $hB_chr);
	}
	if ($uniquexb){
		$refarray_chrB	=	compare_chr_distrib($hB_chr, $hA_chr);
	}
}

# we collect all chromosome-specific transcript in a hash
my %h_diffchr=();


# Process 2 hashes
foreach my $chrA (keys %{$hA_chr} ) {

	# for common, we collect response in array @inter of hashes such as
	# 	$VAR1 = [
	#           {
	#             'ENSCAFT00000000121' => 'ENSCAFT00000000121',
	#             'ENSCAFT00000000971' => 'ENSCAFT00000000971',
	# 		}
	# 		]			
	my @inter;
	
	# for unique of A and/or B
	# just a hash of tx ids
	my %h;
	
	# Test if chrA not in B so we put that in %h_diffchr
	if (!$common && $uniquexa && grep $_ eq $chrA, @{$refarray_chrA} ){
		$h_diffchr{$chrA} = $hA_chr->{$chrA};
		next;
	}
	
	# 2n chromosome
	foreach my $chrB (keys %{$hB_chr} ) {

		# Test if chrB not in A so we put that in %h_diffchr
		if (!$common && $uniquexb && grep $_ eq $chrB, @{$refarray_chrB} ){
			$h_diffchr{$chrB} = $hB_chr->{$chrB};
			next;
		}
		
		next if ($chrA ne $chrB);
                
		# start fork
	    my $pid = $pm->start($chrB) and next;
     	
     	# get ref on hash per chromosome
     	my $refh1 = $hA_chr->{$chrA};
        my $refh2 = $hB_chr->{$chrB};        
        
        
        # Launch interactionCollection
		print STDERR "$chrB\n" if ($verbosity > 0);
		# Launch Intersect
		if (!$common){
			%h = Intersect::Intersect2HsplitFork_compare($refh1, $refh2, $strandedmode, $uniquexa, $uniquexb, $common, $fraction, $verbosity);
		} else {
			@inter = Intersect::Intersect2HsplitFork_compare($refh1, $refh2, $strandedmode, $uniquexa, $uniquexb, $common, $fraction, $verbosity);
		}
    }
    # end fork
    if (!$common){
		$pm->finish(0, \%h);
	} else {
	    $pm->finish(0, \@inter);
	}
}

# wait all sub process
$pm->wait_all_children;


print STDERR "\n";
# Merge 2 hashes of file-specific transcripts if !$common
my %final_diff;
if (!$common){
	%final_diff  = (%retrieved_responses, %h_diffchr);
} 


########
# OUTPUT
#########
# If option $common, we print a double GTF
if ($common){
	# Parse good transcripts by chr
	foreach my $chr (keys %retrieved_responses ) {
		# in option $common $retrieved_responses{$chr} is a array ref
		ExtractFromHash::printDoubleGTF($hA_chr->{$chr}, $hB_chr->{$chr}, $retrieved_responses{$chr}, 'all',  $verbosity)
	}
} else { # else we print simple GTF
	foreach my $chr (keys %final_diff ) {
		ExtractFromHash::printGTF($final_diff{$chr}, 'all',  $verbosity)
	}

}





# from 2 hash ref of split annotation => return an array which contains the list of A chrom that are not in B
sub compare_chr_distrib{
	my ($refha, $refhb) = @_;
	
    my %cmp = map { $_ => 1 } keys %{$refha};
    
	for my $key (keys %{$refhb}) {
        delete $cmp{$key} if exists $cmp{$key};
	}
	my @arkey_justinA = keys %cmp;

    return \@arkey_justinA;
}

__END__

=head1 NAME

intersectRetouch.pl : compare 2 gtf files for matching transcripts (same strand, exon overlap)  or file-specific transcripts...

=head1 SYNOPSIS

intersectRetouch.pl  -a infileA.gtf -b infileB.gtf [options]

=head1 OPTIONS

=over 8

=item B<-a|afile>

Input A .gtf file. [mandatory]

=item B<-b|bfile>

Input B .gtf file. [mandatory]

=item B<-c|common>

write matching transcript(s) between  A and B files... [default TRUE ]

- restricted by -f

- print a double GTF in output

=item B<-ea|uniquexa>

write transcript(s) from A NOT matching transcripts from B ... [default FALSE ]

- restricted by -c

- use --no-common to apply this option

=item B<-eb|uniquexb>

write transcript(s) from B NOT matching transcripts from A ... [default FALSE ]

- restricted by -c

- use --no-common to apply this option

=item B<-f|fraction>

Minimum overlap required as a fraction of A [default : 1 (100% overlap) ].

=item B<-s|strandmode>

Strand mode of overlap: should be 0 (unstranded mode comparison), 1 (same strand) or -1 (for opposite strand e.g antisense) [default : 0 (unstranded) ].

=item B<-p|proc>

number of process [default 1].

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
