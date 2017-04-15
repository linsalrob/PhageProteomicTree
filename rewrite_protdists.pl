#!/usr/bin/perl

# some protdist files don't have all the results on one line anymore, they are split across many lines
# we concatenate all lines

use strict;
my $dir=shift || die "$0 <protdist dir> <dir to write to>";
my $odir = shift || die "$0 <protdist dir> <dir to write to>";
mkdir $odir, 0755 unless (-e $odir);

opendir(DIR, $dir) || die "can't open $dir";
foreach my $file (grep {$_ !~ /^\./} readdir(DIR)) {
	open(IN, "$dir/$file") || die "can't open $file for reading";	
	open(OUT, ">$odir/$file") || die "can't open $file for writing";
	my $first;
	while (<IN>)
	{
		chomp;
		unless ($first) {
			print OUT;
			$first=1;
			next;
		}
		if (/^\d/) {print OUT "\n$_"}
		else {print OUT}
	}
	print OUT "\n";
	close IN;
	close OUT;
}

	
