#!/usr/bin/perl -w

use strict;
use warnings;
use Utils;
use Data::Dumper;
use Tie::IxHash;
use Parser;


my $infile=shift;

open INFILE, $infile or die "cannot open $infile\n";


my %h_gene = Parser::parseGtf($infile,10);

# foreach my $line  (<INFILE>){
# 	chomp $line;
# 	if ($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)"; (\w.*)$/ || 	
# 	$line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+gene_id "(\S+)"; transcript_id "(\S+)";$/){
# 		my $chr     	=   $1; 
# 		my $biotype		=   $2; 
# 		my $level		=   $3; 
# 		my $beg_feat    =   $4; 
# 		my $end_feat    =   $5;  
# 		my $strand      =   $7; 
# 		my $gene_id		=   $9;
# 
# 		$h_gene{$gene_id}->{"chr"}			=   $chr;
# 		$h_gene{$gene_id}->{"startg"}		=   Utils::min2($h_gene{$gene_id}->{"startg"}, $beg_feat);
# 		$h_gene{$gene_id}->{"endg"}		=   Utils::max2($h_gene{$gene_id}->{"endg"}, $end_feat);	
# 	}
# }
print Dumper \%h_gene;


tie my %h_gene_sorted, 'Tie::IxHash';

foreach my $key (sort { $h_gene{$a}->{'startt'} <=> $h_gene{$b}->{'startt'} }  keys %h_gene){

	my $value = $h_gene{$key};
    printf( "%s\t%s\t%s\n", $key, $value->{'chr'}, $value->{'startt'});
	$h_gene_sorted{$key}=$h_gene{$key};
}

# define the ordered hash
tie my %h_transcript_ordered, 'Tie::IxHash';

my %h_transcript_ordered = Parser::orderHashStart(\%h_gene, 0);
print Dumper \%h_transcript_ordered;

die;




print Dumper \%h_gene_sorted;


close INFILE;