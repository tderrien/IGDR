#!/usr/bin/perl -w

#
# May 2014
#
# tderrien@univ-rennes1.fr
# adapted from script here:
# 	https://www.biostars.org/p/2226/
# Split a multi fasta file
#
# run 'split_fasta.pl -h' for help
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
# output prefix file
my $prefix 	= 	"split";
# number of sequences per output file(s)
my $number	= 	1 ;
# help and man
my $help	=	0;
my $man		=	0;


################################ GetOptions###############################################
GetOptions (
	"infile|i=s" 	=>	\$infile,
	"prefix|p=s"	=>	\$prefix,
	"number|n=i"	=>	\$number,
	"help|?"	=>	\$help,
	"man|m"		=>	\$man
) or pod2usage(2);	

# Print help if needed
pod2usage(1) if $help;
pod2usage(-exitval=>0, -verbose=>2) if $man;

# Mandatory options
pod2usage ("-i option: I need a non-empty input file...\n") unless (-s $infile);


my $in  = new Bio::SeqIO(-file  => $infile);

my $count = 0;
my $fcount = 0;
my $out ;
while (my $seq = $in->next_seq) {
    if ($count % $number == 0) {
        $fcount++;
        $out = new Bio::SeqIO(-file => ">$prefix.$fcount", -format=>'fasta');
    }
    $out->write_seq($seq);
    $count++;
}


__END__

=head1 NAME

split_fasta.pl - Split a multi-fasta file in files of n sequence(s)

=head1 SYNOPSIS

perl split_fasta.pl -i <MultiFasta_File> [Options]

	Required parameters
		-i <MultiFasta_File>
	
	Optional Parameters
		-p prefix of the output files			[default 'split']
		-n number of sequence per file 			[default '1']
		-h help message					[default '0']
		-m man message					[default '0']

=head1 Optional parameters

=item B<--prefix (or -p)>

The string prefix of the output files (default 'split')

=item B<--number (or -n)>

The desired number of sequence per file file(s) (default '1')

=item B<--help (or -h)>

This Help message

=item B<--man (or -m)>

Th manual page

=cut

