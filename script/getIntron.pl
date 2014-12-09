#!/usr/bin/perl -w

use strict;
use warnings;

use Data::Dumper;
use Parser;
use ExtractFromHash;

my $infile=shift;

my %h = Parser::parseLevelGtfhashKey($infile);
my %hi = ExtractFromHash::getIntronFromGtfHash(\%h);
ExtractFromHash::printGtfHash(\%h);

