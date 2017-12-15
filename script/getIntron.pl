#!/usr/bin/perl -w


# TD:
# AIM: get introns from a .gtf file having exon informations
#############################################################

use strict;
use warnings;


use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Parser;

use Getopt::Long;
use Pod::Usage;
use File::Basename;
use Data::Dumper;

# Variables
my $progname=basename($0);

# Variables
my $infile		='';
my %filtertag;
my $man 		= 0;
my $help 		= 0;
my $verbosity	= 0;


## or if usage was explicitly requested.
GetOptions(
	'i|infile=s'	 	=> \$infile,
	"f|filter=s"		=> \%filtertag,
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'man' 				=> \$man
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test parameters
pod2usage("Error: Cannot read your input GTF file '$infile'...\nFor help, see:\n$progname --help\n") unless( -r $infile);


my $h				= Parser::parseGTF($infile, 'exon',  undef , \%filtertag , $verbosity);
my $refintrons		= undef;

foreach my $tr ( keys %$h ) {

    $refintrons                 =        ExtractFromFeature::getIntrons($h->{$tr}->{'feature'});

	foreach my $feat (@$refintrons) {


		print join("\t",$h->{$tr}->{"chr"}, $h->{$tr}->{"source"}, $feat->{"feat_level"}, $feat->{"start"}, $feat->{"end"}, $h->{$tr}->{"score"}, $h->{$tr}->{"strand"}, ".");
		print "\tgene_id \"".$h->{$tr}->{"gene_id"}."\"; transcript_id \"$tr\";\n";

	}
}




__END__

=pod

=encoding UTF-8

=head1 NAME

getIntron.pl.pl - get introns from a .gtf file having exons levels

=head1 VERSION

version 0.01

=head1 SYNOPSIS

getIntron.pl.pl -i INPUT.gtf 

=head1 DESCRIPTION

get introns of a .GTF file w.r.t to different user's filter

=head1 OPTIONS

=head2 General

  --help                Print this help
  --man                 Open man page
  --verbosity		Level of verbosity
  

=head2 Mandatory arguments

  -i,--infile=file.gtf		Specify the GTF file in input

=head2 Filtering arguments

  -f,--filter			filtering attributes such as key=val (-f chr=X,34 -f transcript_biotype=lincRNA,antisense) [default: undef]
  
  	

=head1 AUTHORS

=over 4

=item -
Thomas DERRIEN <tderrien@univ-rennes1.fr>


=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2015 by IGDR - CNRS

=cut
