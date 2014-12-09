#!/usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Parser;

use Bio::EnsEMBL::Registry;
use Data::Dumper;
use EnsemblAPI;


# input  file 
my $vcffile='';
my $gtffile='';

# Connecxion to emsembl API
my $reg			=	"Bio::EnsEMBL::Registry";
# target species 
my $target		=	"";
# orthology type
#my $orthtype	=	"one2one";
# verbosity, help and man
my $verbosity	=	0;
my $help		=	0;
my $man			=	0;

    
# Parsing parameters
my $res=GetOptions (
    'vcf|v=s'           =>	\$vcffile,
	'gtffile|g=s'		=>	\$gtffile,
	'verbosity' 		=> \$verbosity,
	'help|?'			=> \$help,
	'man|m'				=>	\$man,
) or pod2usage ("Try '$0 --help' for more information");


# Parsing parameters 
if ($verbosity>0){ print STDERR "Parsing parameters...\n";}

pod2usage ( -verbose => 1) if $help;
pod2usage ( -verbose => 2) if $man;
pod2usage ("I need an input .vcf file (option --vcf or -v)...\n") unless (-r $vcffile);
pod2usage ("I need an input .gtf file (option -g)...\n") unless (-r $gtffile);



print STDERR "Checking options end ...\n";

# Intersect command
my $command="closestBed -a $vcffile -b $gtffile -d -t first |";
print STDERR "$command\n";
open(IN, $command) || die "can't run $command";

my %h =();
my %gene_iddup = ();

while (<IN>){
    
    next if(/^##/); #ignore header
    chomp;
    
    #initialize hash with possible missing keys
    my %attribs = ('gene_name' => 'NA', 'gene_biotype' => 'NA');

    # parse line
    my ($chr, $pos, $snpid, $ref, $alt, $qual, $filter, $info, 
    $chrgtf, $source, $type, $start, $end, $score, $strand, $phase, $attributes, $distance) = split("\t");

    # store SNP data
    $h{$snpid}->{"chr"}			=   $chr;
	$h{$snpid}->{"pos"}	    	=   $pos;
	$h{$snpid}->{"REF"}	    	=   $ref;
	$h{$snpid}->{"ALT"}         =   $alt;
			

    # last element of an array
    $attribs{'distance'} = $distance;
    
    
    # parse gtf attributes
    my @add_attributes = split(";", $attributes);
    
    # store ids and additional information in second hash
    foreach my $attr ( @add_attributes ) {
        next unless $attr =~ /^\s*(.+)\s(.+)$/;
        my $c_type  = $1;
        my $c_value = $2;
        $c_value=~ s/\"//g;
        
        if($c_type  && $c_value){
            $attribs{$c_type} = $c_value;
        }

    }
    push (@{$h{$snpid}->{"OVERLAP"}}, \%attribs);     # if (%attribs) tests whether %attrib is not empty
}

# Extract unique gene data overlap in hash '%genes'
#   - gene_id
#   - gene_biotype
#   - gene_name
#     print $snp, scalar @{$h{$snp}->{"OVERLAP"}}, "\n";
#######################################
# print Output
for my $snp (keys %h){

    print $h{$snp}->{"chr"},"\t",
	$h{$snp}->{"pos"},"\t",
	$h{$snp}->{"REF"},"\t",
    $h{$snp}->{"ALT"},"\t",
	$snp,"\t";
    
    foreach my $feat (@{$h{$snp}->{"OVERLAP"}}) {
    
 	    print " ",$feat->{"gene_id"},";",$feat->{"gene_name"},";",$feat->{"gene_biotype"},";",$feat->{"distance"};
    }
    print "\n";
}    
# 
