#!/usr/bin/perl -w

#
# Jan 2015
#
# tderrien@univ-rennes1.fr
# get length/size of a .mfa sequence
# adapted from http://www.bioperl.org/wiki/HOWTO:SeqIO
# run 'prog.pl -h' for help
#
##########################################################################################


############################ library #####################################################
use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

############################# Parameters #################################################


# mfasta file
my $infile	= "";

# help and man
my $help	=	0;
my $man		=	0;


################################ GetOptions###############################################
GetOptions (
	"infile|i=s" 	=>	\$infile,
	"help|?"	=>	\$help,
	"man|m"		=>	\$man
) or pod2usage(2);	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-exitval=>0, -verbose=>2) if $man;

# Mandatory options
pod2usage ("-i option: I need a non-empty input file...\n") unless (-s $infile);


my $in  = new Bio::SeqIO(-file  => $infile);
my @seq_array;

while( my $seq = $in->next_seq() ) {
    push(@seq_array,$seq);
}
 
# now do something with these. First sort by length,
# find the average and median lengths and print them out
@seq_array = sort { $a->length <=> $b->length } @seq_array;
 
my $total = 0;
my $count = 0;
foreach my $seq ( @seq_array ) {
    print $seq->id,"\t",$seq->length,"\n";
    $total += $seq->length;
    $count++;
}
 
print STDERR "# Mean length ",int($total/$count)," nt Median ",$seq_array[$count/2]->length," nt\n";


__END__

=head1 NAME

fasta_size.pl - print size of fasta sequences (and mean and mediane size in STDERR)

=head1 SYNOPSIS

perl fasta_size.pl -i <MultiFasta_File> [Options]

	Required parameters
		-i <MultiFasta_File>
	
	Optional Parameters

		-h help message					[default '0']
		-m man message					[default '0']

=head1 Optional parameters


=item B<--help (or -h)>

This Help message

=item B<--man (or -m)>

Th manual page

=cut

