#!/usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Parser;

use Bio::EnsEMBL::Registry;
use Data::Dumper;
use EnsemblAPI;

#########################
# TODO : get human genes by API 
#		: Need API to be working on ensembl_compara in genocluster mysql.
#########################
# my $ensg = EnsemblAPI::getOrtholog('ENSG00000012048','canis_familiaris');
# print "$ensg\n";
# die;

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

# my $length='';

# GetOptions ('length|height=f' => \$length) or pod2usage ("Try '$0 --help' for more information");
    
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
#pod2usage ("I need a target species e.g Mus_musculus...\n") unless ($target ne "");
#pod2usage ("I need a valid orthology type...\n") unless ($orthtype eq "one2one" || $orthtype eq "one2many" ||$orthtype eq "all");

# my %htest=Parser::parseGtf($gtffile);
# print Dumper %htest;



print STDERR "Checking options end ...\n";

# Intersect command
my $command="intersectBed -a $vcffile -b $gtffile -wa -wb -loj |";
print STDERR "$command\n";
open(IN, $command) || die "can't run $command";

my %h =();
my %gene_iddup = ();

while (<IN>){
    
    next if(/^##/); #ignore header
    chomp;
    
    my %attribs = ('gene_name' => 'NA', 'gene_biotype' => 'NA', 'gene_id' => 'NA');

    # parse line
    my ($chr, $pos, $snpid, $ref, $alt, $qual, $filter, $info, 
    $chrgtf, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");

    # store SNP data
    $h{$snpid}->{"chr"}			=   $chr;
	$h{$snpid}->{"pos"}	    	=   $pos;
	$h{$snpid}->{"REF"}	    	=   $ref;
	$h{$snpid}->{"ALT"}         =   $alt;
			
    # parse gtf attributes
    my @add_attributes = split(";", $attributes);
    # store ids and additional information in second hash
    foreach my $attr ( @add_attributes ) {
        next unless $attr =~ /^\s*(.+)\s(.+)$/;
        my $c_type  = $1;
        my $c_value = $2;
        $c_value=~ s/\"//g;
        
#         if(exists($attribs{$c_type})){
#             warn "$c_type key already exists with ",$attribs{$c_type}," value...\n";
#         }
        if($c_type  && $c_value){
            $attribs{$c_type} = $c_value;
        }
    }

    # count how many time the pair *snp_id.gene* (i.e $attribs{'gene_id'}) is found
    $gene_iddup{$snpid.$attribs{'gene_id'}}++;

    # add the ref to hash only for unique key *snp_id.gene* in order to remove duplicates
    if ($gene_iddup{$snpid.$attribs{'gene_id'}}==1){
        push (@{$h{$snpid}->{"OVERLAP"}}, \%attribs);     # if (%attribs) tests whether %attrib is not empty
    }
}

# print Dumper \%h;

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
    
 	    print " ",$feat->{"gene_id"},";",$feat->{"gene_name"},";",$feat->{"gene_biotype"};
    }
    print "\n";
}    

