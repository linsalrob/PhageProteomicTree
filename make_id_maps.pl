#!/usr/bin/perl -w

use strict;

foreach my $f (qw[genome_id.map protein_genome_lengths.txt]) {
	die "$f exists, not overwriting" if (-e $f);
}

my ($idmap, $real, $gennamef)=@ARGV;
unless ($idmap && $real && $gennamef) {	
	die "Usage: $0 <id.map> <tuple of real id, genome id, length> <tuple of genomeid, real name, abbreviated name for tree>";
}




my %id;
open(IN, $idmap) || die "Can't open $idmap";
while (<IN>) {
	chomp;
	my @a=split /\t/;
	$id{$a[0]}=$a[1];
}
close IN;

my %genomename;
open(IN, $gennamef) || die "can't open $gennamef";
while (<IN>) {
	chomp;
	my ($id, $full, $abbr)=split /\t/;
	if (defined $abbr) {
		$genomename{$id}=$abbr;
	} else {
		print STDERR "Warning. No abbreviation for $id, so using $id as the name!\n";
		$genomename{$id}=$id;
	}
}

open(IN, $real) || die "Can't open $real";
open(OUT, ">protein_genome_lengths.txt") || die "Can't write to protein_genome_lengths.txt";
my %genome; my $g=0;
while (<IN>) {
	chomp;
	my ($prot, $gen, $len)=split /\t/;
	unless (defined $genome{$gen}) {$genome{$gen}=++$g}
	unless ($id{$prot}) {die "no id for $prot"}
	print OUT join("\t", $genome{$gen}, $id{$prot}, $len), "\n";
}
close OUT;
	
	
	
open(GEN, ">genome_id.map") || die "Can't open genome_id.map";
foreach my $gen (sort {uc($genomename{$a}) cmp uc($genomename{$b})} keys %genome) {
	unless (defined $genome{$gen}) {
		print STDERR "No genome for $gen ($genomename{$gen})\n";
		$genome{$gen} = "XXX";
	}
	unless (defined $genomename{$gen}) {
		                print STDERR "No genome for $gen ($genome{$gen}\n";
				$genomename{$gen} = "XXX";
			}
	print GEN $genome{$gen}, "\t", $genomename{$gen}, "\n";
}
close GEN;
