#!/usr/bin/perl -w

use strict;
my $genomenamef=shift || die "genome name file?";
my $treefile=shift || die "tree file?";
my %gn;
open(IN,  $genomenamef) || die "can't open  $genomenamef";
while (<IN>)
{
	chomp;
	my @a=split /\t/;
	$a[1] =~ s/\s+/_/g;
	$a[1] =~ s/\W//g;
	$gn{"gnm".$a[0]}=$a[1];
}
close IN;

open(IN, $treefile) || die "can't open $treefile";
while (<IN>)
{
	while (/(gnm\d+)\D/)
	{
		my $h = $1;
		my $m = $gn{$h};
		if ($m) {
			s/$h(\D)/$m$1/g;
		}
		else {
			print STDERR "No match for $1\n";
			my $H=uc($h);
			s/$h(\D)/$H$1/g;
		}
	}

	print;
}
close IN;

