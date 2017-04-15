#!/usr/bin/perl -w

# convert proteins.faa to a list of protein IDs, genome IDs, and protein length
#

use strict;
use Rob;
my $rob = new Rob;

my $file = shift || die "fasta file of proteins. Probably proteins.faa";
my $fa=Rob->read_fasta($file);
map {
	m/^(\S+)\s+\[(.*?)\]/;
	if ($1 && $2) {print join("\t", $1, $2, length($fa->{$_})), "\n"}
	else {print STDERR "Can't parse $_\n"}
} sort keys %$fa;

