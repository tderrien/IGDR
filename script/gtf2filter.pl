#!/usr/bin/perl -w



# Perl libs
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Basename;


# Own lib : /home/genouest/umr6061/recomgen/tderrien/bin/ThomasPerl/lib/
use Parser;
use ExtractFromHash;
use ExtractFromFeature;
use Intersect;
use Utils;



my $progname=basename($0);

# Variables
my $infile		='';
my %filtertag;
my $man 		= 0;
my $help 		= 0;
my $verbosity	= 0;

my $longest		=	0; # extract longest tx per locus

my $minsize         = 200;
my $biexonicsize	= 25;
my $monoexonic      = 0; # -1 keep monoexonicAS, 1 keep all monoexonic, 0 remove all monoexonic
						 # restricted by $linconly

## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
	'i|infile=s'	 	=> \$infile,
    's|size=i' 			=> \$minsize,
    'biex=i' 			=> \$biexonicsize,    
    'monoex=i' 			=> \$monoexonic,
	"f|filter=s"		=> \%filtertag,
	'longest!'			=> \$longest,
	'v|verbosity=i'		=> \$verbosity,
	'help|?' 			=> \$help,
	'man' 				=> \$man
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;


# Test parameters
pod2usage("Error: Cannot read your input GTF file '$infile'...\nFor help, see:\n$progname --help\n") unless( -r $infile);
pod2usage ("- Error: --monoex option '$monoexonic' should be 1 keep all monoexonic, 0 remove all monoexonic \n") if ($monoexonic != 0 and $monoexonic != 1 );

#############################################################


# Parsing candidate lncRNAs
my $reflnc		= Parser::parseGTF($infile, 'exon',  undef , \%filtertag , $verbosity);

 

my ($ctminsize, $ctmonoexonic, $ctdubious) = (0,0,0);
# Filtering steps
foreach my $tx (keys %{$reflnc}){

    my $txfeatures  =   $reflnc->{$tx}->{'feature'};
    
    # size
    my $size = ExtractFromFeature::features2size($txfeatures, 0);
    if ($size <$minsize){
        $ctminsize++;
        delete $reflnc->{$tx};
    }
    
    # nb exon
    my $nbexon = ExtractFromFeature::features2nbExon($txfeatures, 0);
    if (!$monoexonic && $nbexon == 1){
        $ctmonoexonic++;
        delete $reflnc->{$tx};    
    }
    
    # dubious biexonic
    my $dubious = ExtractFromFeature::features2dubExon($txfeatures, $biexonicsize, 0);
    if ($dubious && $nbexon == 2){
        $ctdubious++;
        delete $reflnc->{$tx};     
    }
    
    
}

print STDERR "> Filter size ($minsize): $ctminsize\n";
print STDERR "> Filter monoexonic ($monoexonic): $ctmonoexonic\n";
print STDERR "> Filter biexonicsize ($biexonicsize): $ctdubious\n";
print STDERR ">> Transcripts left after fitler(s): ",scalar keys (%{$reflnc}),"\n";

# Get longest transcripts per gene
if ($longest) {
	print STDERR "Extract longest tx per locus : $longest\n";
	my $reflnc_long	=	ExtractFromHash::getCumulSizeFromGtfHash($reflnc,  $verbosity, $longest);

	ExtractFromHash::printGTF($reflnc_long, 'all',  $verbosity);
	
} else {

	ExtractFromHash::printGTF($reflnc, 'all',  $verbosity);

}




__END__

=pod

=encoding UTF-8

=head1 NAME

gtf2filter.pl - Extract, filter gtf file

=head1 VERSION

version 0.01

=head1 SYNOPSIS

gtf2filter.pl -i INPUT.gtf 

=head1 DESCRIPTION

Filter a .GTF file w.r.t to different user criteria 

=head1 OPTIONS

=head2 General

  --help                Print this help
  --man                 Open man page
  --verbosity		Level of verbosity
  

=head2 Mandatory arguments

  -i,--infile=file.gtf		Specify the GTF file to be filtered 
  

=head2 Filtering arguments

  -f,--filter			filtering attributes such as key=val (-f chr=X,34 -f transcript_biotype=lincRNA,antisense) [default: undef]
  -s,--size=200			Keep transcript with a minimal size (default 200)
  --monoex=0|1			Keep monoexonic transcript(s): mode to be selected from :  1 keep all monoexonic, 0 remove all monoexonic	[default 0]
  --biex=25			Discard biexonic transcripts having one exon size lower to this value (default 25)
  --longest			get longest tx per locus (default 'FALSE')
  
  	

=head1 AUTHORS

=over 4

=item -
Thomas DERRIEN <tderrien@univ-rennes1.fr>


=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2014 by IGDR - CNRS

=cut
