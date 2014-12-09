#!/usr/bin/perl 

use strict;
use warnings;
use Parser;
use ExtractFromHash;
use Data::Dumper;


my $infile=shift;

my $refh= Parser::parseGTF($infile,'exon', 1);

print Dumper $refh;
 ExtractFromHash::printGtfHash($refh, 'all', 10);


