#!/usr/bin/perl -w



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

# Extract 10 first sequences by default
my $rangefilter 	=	undef;
my ($min, $max)= 0,0;



my $minsize         = 0;
my $biexonicsize	= 0;
my $monoexonic      = 1; # -1 keep monoexonicAS, 1 keep all monoexonic, 0 remove all monoexonic
						 # restricted by $linconly

#my $has_complete_ORF	=	0 # if 1, will only get transcript having start and stop codon annotated (some transcript do not have cds_{start,end}_NF while they do not start or finish by start or stop e.g ENST00000429039 ih GRCH37/hg19)
								# see /home/genouest/umr6061/recomgen/dog/script/FEELnc++/data/human/randomTxperGene/without_cds_end-start_NF/error.id

my $longest		=	0; # extract longest tx per locus
my $random		=	0; # extract one random tx per locus


## Parse options and print usage if there is a syntax error,
## or if usage was explicitly requested.
GetOptions(
	'i|infile=s'	 	=> \$infile,
    's|size=i' 			=> \$minsize,
    'biex=i' 			=> \$biexonicsize,    
    'monoex=i' 			=> \$monoexonic,
	"f|filter=s"		=> \%filtertag,
	"n|number=s"		=> \$rangefilter,
	'longest!'			=> \$longest,
	'r|random!'			=> \$random,
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

my $filt_reflnc;
if (defined $rangefilter){

	if (index($rangefilter, '-') != -1) { # if range
		($min, $max) = split /-/, $rangefilter;
		if ($max<$min){my $temp=$max; $max=$min;$min=$temp}
		print STDERR "Extracting txs by NUMBER from: $min->$max\n";

		# get tx_ids from hash in array
		my @k = keys(%{$reflnc});
		# subselect from min to max tx
		my @sk= @k[$min .. $max];
		%{$reflnc} = map { $_ => $reflnc->{$_}} @sk;

	}
}





# Get longest transcripts per gene
if ($longest) {
	print STDERR "Extract longest tx per locus : $longest\n";
	my $reflnc_long	=	ExtractFromHash::getCumulSizeFromGtfHash($reflnc,  $verbosity, $longest);

	ExtractFromHash::printGTF($reflnc_long, 'all',  $verbosity);
	
} elsif ($random) {

	&printGTFrdmTx($reflnc, 'all',  $verbosity);

    
} else {

	ExtractFromHash::printGTF($reflnc, 'all',  $verbosity);

}



# only print one random transcript per gene
sub printGTFrdmTx{
	# parsing a hash in sub need dereference in shift
	my ($refh, $mode, $verbosity)	= @_;
	$mode 			||= "all"; 		# number of fields to be printed
	$verbosity 		||= 0;

	# print if verbosity	
	print STDERR "\nPrinting transcripts GTF...\n" if ($verbosity > 0);
	
	my %seen_gn;
	
# 	print Dumper $refh;
	# Parse gtfHash to be printed
	
	for my $tr (keys %{$refh}){

        # increment hash gene counter
        $seen_gn{$refh->{$tr}->{"gene_id"}}++;
        
        # next if already seen gene
        next if ($seen_gn{$refh->{$tr}->{"gene_id"}} > 1 );
	    
	    
	    # print normally as printGTF
		foreach my $feat1 (@{$refh->{$tr}->{"feature"}}) {
			print join("\t",$refh->{$tr}->{"chr"}, $refh->{$tr}->{"source"}, $feat1->{"feat_level"}, $feat1->{"start"}, $feat1->{"end"}, $refh->{$tr}->{"score"}, $refh->{$tr}->{"strand"}, $feat1->{"frame"});
 			print "\tgene_id \"".$refh->{$tr}->{"gene_id"}."\"; transcript_id \"$tr\";";
 			
			if ($mode eq "all"){
				my %tmph = %{$feat1};
				# delete unnecesserary keys from hash %tmph
				delete @tmph{qw/feat_level start end strand frame/};
				for (sort keys %tmph){
					print " $_ \"$tmph{$_}\";" if (defined $tmph{$_}); # id defined in order to test if parsinf $extrafield is OK (i.e in case people select FPKM and it does not exists in the file)
				}
			}	
			print "\n";
		}		
	}
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
  -n,--number			range of tx to be extract (-n 0-9 extract 10 first transcripts... WARNINGs with hash order) [default undef]
  -s,--size=0			Keep transcript with a minimal size (default >0)
  --monoex=0|1			Keep monoexonic transcript(s): mode to be selected from :  1 keep all monoexonic, 0 remove all monoexonic	[default 1]
  --biex=0			Discard biexonic transcripts having one exon size lower to this value (default >0)
  --longest			get longest tx per locus ... following filters (default 'FALSE')
  -r,--random 			get random tx per locus ... following filters (default 'FALSE')
  

$ perl extract_fasta.pl -i infile.fa  B<-f 1-10>

: will extract 10 first sequences from infile.fa  	

=head1 AUTHORS

=over 4

=item -
Thomas DERRIEN <tderrien@univ-rennes1.fr>


=back

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2014 by IGDR - CNRS

=cut
