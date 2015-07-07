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
my $minfrac_over 	= 0 ;
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
	'va|uniquexa!' 		=> \$uniquexa,
	'f|fraction=f' 		=> \$minfrac_over,	
	's|strandmode=i' 	=> \$strandedmode,
	'p|proc=i' 			=> \$proc,
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'man' 				=> \$man
) or pod2usage(2);	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-exitval=>0, -verbose=>2) if $man;

print STDERR "Options:\n- infileA: $infileA\n- infileB: $infileB\n- uniquexa: $uniquexa\n- minfrac_over: $minfrac_over\n- strandedmode: $strandedmode\n" if ($verbosity > 5);

# Test parameters
pod2usage("Error: Cannot read your input A .gtf file '$infileA'...\nFor help, see:\n$progname --help\n") unless( -r $infileA);
pod2usage("Error: Cannot read your input B .gtf file '$infileB'...\nFor help, see:\n$progname --help\n") unless( -r $infileB);
pod2usage ("- Error: \$minfrac_over option '$minfrac_over' should be a float between 0 and 1 [0-1] (e.g 0.5 if 50% overlap is required)\n") unless ($minfrac_over >= 0 and $minfrac_over <= 1);

# pod2usage("Error: Cannot report common tx while unique A or B activated. Use '--no-common' option to exclude matching tx.\n\nFor help, see:\n$progname --help\n") if( $common && $uniquexa);
pod2usage ("- Error: \$strandedmode option '$strandedmode' should be 0 (unstranded mode comparison), 1 (same strand) or -1 (for opposite strand e.g antisense)\n") if ($strandedmode !=0 && $strandedmode !=1 && $strandedmode != -1 );


my $splitbychr	=	1;
my $parseextraf	=	undef;
my %biotype;
# Parsing
my $hA_chr		= Parser::parseGTF($infileA, 'exon',  $splitbychr, \%biotype , $verbosity);
my $hB_chr		= Parser::parseGTF($infileB, 'exon',  $splitbychr, \%biotype , $verbosity);


# Launch parralel
my $pm = Parallel::ForkManager->new($proc);
my %refhchild = ();  # for collecting responses
$pm -> run_on_finish (
  sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $result_ref) = @_;
		my $chr = $ident;

		warn("Child $chr killed by signal $exit_signal"), return if $exit_signal;
		warn("Child $chr exited with error $exit_code"),  return if $exit_code;
		warn("Child $chr encountered an unknown error"),  return if !$result_ref;

	   $refhchild{$chr} = $result_ref;
  }
);

# For option uniqueA or B, we need to check whether some txs are specific of ONE chromosome
# that will be miss by the chr to chr comparison
my $refarray_chrA;
if ($uniquexa){
	$refarray_chrA	=	compare_chr_distrib($hA_chr, $hB_chr);
}

# we collect all chromosome-specific transcript in a hash
my %h_diffchr=();


# Process 2 hashes
foreach my $chrA (keys %{$hA_chr} ) {

 	# array tmp storing tx IDs to remove, the correct hash is : $refhchild
	# 	$VAR1 = [
	#           {
	#             'ENSCAFT00000000121' => 'ENSCAFT00000000121',
	#             'ENSCAFT00000000971' => 'ENSCAFT00000000971',
	# 		}
	# 		]			
	my @pairs;
	
	# Test if chrA not in B so we put that in %h_diffchr
	if ($uniquexa && grep $_ eq $chrA, @{$refarray_chrA} ){
		$h_diffchr{$chrA} = $hA_chr->{$chrA};
		next;
	}
	
	# 2n chromosome
	foreach my $chrB (keys %{$hB_chr} ) {

		
		next if ($chrA ne $chrB);
                
		# start fork
	    my $pid = $pm->start($chrB) and next; 
        
        # Launch interactionCollection
		print STDERR "$chrB\n" if ($verbosity > 0);
		# Launch Intersect

		# Get all tx from A that intersect with B 
		# 			* wrt to strand and fraction option (i.e restricted by -s)
		# 			* in format such chr => {A1	=> B1}
		#						    	=>  {A2	=> B1}
		# 								=>  .....
		@pairs = Intersect::getIntersect($hA_chr->{$chrA}, $hB_chr->{$chrB}, $strandedmode, $minfrac_over, $verbosity);
    }
    # end fork
	$pm->finish(0, @pairs);
}

# wait all sub process
$pm->wait_all_children;

# print Dumper \%refhchild;

if ($uniquexa){
	# remove matching transcripts from hash of process chr
	foreach my $thread (keys %refhchild){
    	my @arayofhash =  @{$refhchild{$thread}};
#     	print Dumper \@arayofhash;
    	my @idstorm;
    	for my $ref (@arayofhash){
			push @idstorm, keys %{$ref};
    	}
# 		print STDERR "Filter because  $_ - $idstorm[$_]...\n" for (@idstorm);
    	
    	delete @{$hA_chr->{$thread}}{@idstorm} ;  
    }
}


# Merge 2 hashes of file-specific transcripts if !$common
my %final_diff;
if ($uniquexa){
	%final_diff  = (%{$hA_chr}, %h_diffchr);
} 



# ########
# # OUTPUT
# #########
########
# OUTPUT
#########
# If option $common, we print a double GTF
if (!$uniquexa){
	# Parse good transcripts by chr
	foreach my $chr (keys %refhchild ) {
		# in option $common $retrieved_responses{$chr} is a array ref
		ExtractFromHash::printDoubleGTF($hA_chr->{$chr}, $hB_chr->{$chr}, $refhchild{$chr}, 'all',  $verbosity)
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


- restricted by -c

- use --no-common to apply this option

=item B<-f|fraction>

Minimum overlap required as a fraction of A [default : 0 (1 nt overlap) ].

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
