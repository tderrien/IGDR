#!/bin/perl -w

use strict;
use warnings;
use Parser;
use Data::Dumper;

my $infile = shift;


my %h_chr = Parser::parseGTFSplitchr($infile,10);

print Dumper \%h_chr;


