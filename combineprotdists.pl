#!/usr/bin/perl -w

#    Copyright 2001, 20002 Rob Edwards
#    For updates, more information, or to discuss the scripts
#    please contact Rob Edwards at redwards@utmem.edu or via http://www.salmonella.org/
#
#    This file is part of The Phage Proteome Scripts developed by Rob Edwards.
#
#    Tnese scripts are free software; you can redistribute and/or modify
#    them under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    They are distributed in the hope that they will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    in the file (COPYING) along with these scripts; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#



# combinprotdists.pl

# a trimmed down version that doesn't complain if the old data and new data don't match. it works fine, but uses less memory

# This version doesn't rely on mysql - you need to provide some text files for the mapping.

# new version. This will assign a distance of $penalty to any sequence that does not match, and to all empty
# spaces. There is a good rationale for this. The The Dayhoff PAM matrix scoring system returns a percent of
# the amino acids that are likely to have changed. Therefore a 100% score means that they have all changed.

# We will make an average, and include the number of scores used to calculate the average. Then we can fitch
# it with the subreplicas option.

# the subreplicate number will be added if the protein appears, but is not similar. (i.e. an average score of 100 2
# means that two proteins were found to be similar but were to distant for a score. But an average score of
# 100 0 means that no proteins were found!

# in this version you can select -n for no padding of missing proteins. Padding runs through the genome and if
# no match to another genome is found it increments the score by 100 for each protein in the query (line) genome
use strict;

my $usage = "combineprotdists.pl <dir of prot dists> <matrix filename> <number of genomes used> <options>\nOPTIONS\n";
$usage .= "\t-d data file : must be in the context [genome, protein id, protein length]\n";
$usage .= "\t-n DON'T pad missing proteins with worst score (supplied with the -p option)\n\t-m print out all protein matches";
$usage .= "\n\t-p # penalty for being bad Rob. Default is 100\n\t-l factor lengths of proteins into the scores\n";
$usage .= "\t-lp penalize based on protein lengths (otherwise it will be whatever the penalty is)\n";
$usage .= "\t-s # skip proteins with match >= this value (treat as if there is no match)\n";
$usage .= "\t-c report progress (will be put on STDERR, but you may want to redirect this\n";
$usage .= "\t-a don't average the protein scores (only works with -l)\n";



my $dir= shift || &niceexit($usage);
my $matrixfilename = shift || &niceexit($usage);
my $nogenomes = shift || &niceexit($usage);

my $args= join (" ", @ARGV);
my $pad=1; my $penalty =100;
my $skip=100000000; # start with skip unreasonably high. Nothing will be bigger than this. Reset if called on command line
my $genomedatafile;
if ($args =~ /-n/) {$pad=0}
if ($args =~ /-p\s+(\d+)/) {$penalty=$1}
if ($args =~ /-s\s+(\d+)/) {$skip=$1; print STDERR "Using skip of $skip\n"}
if ($args =~ /-d\s+(\S+)/) {$genomedatafile=$1}
 &niceexit($usage) unless (defined $genomedatafile);

print STDERR "Using PENALTY of $penalty and PAD of $pad\n";

my %noorfs;
&getnoorfs;

my %proteinlength;
if ($args =~ /-l/) {
my $length = &getprotlengths;
%proteinlength = %$length;
}

my @matrix; 

my @newmatch; my @newmatchcount; my @newmismatch; my @newmismatchcount;

my $genomedata; # the replacement for dbh

my @count; 
my @proteinmatches; my @genematches; my $minmatch = 100; my $maxmatch=1;
my %linegenomecount;
# read each file one at a time, and add the data to an array
opendir(DIR, $dir) || &niceexit("Can't open $dir");
print STDERR "Reading the files\n";
while (my $file=readdir(DIR)) {
   next if ($file =~ /^\./);
   open (IN, "$dir/$file") || &niceexit("Can't open $dir/$file");
   my @genomes = ('0');
   my @genes = ('0');
   my @dists;
   my %genomeshash;
   my %checkdupgenomes;
   # we need to know all the genomes before we can store the data. Therefore
   # read and store each line in @dists
   # then get all the genome numbers and store them in an array
   while (my $line = <IN>) {
      chomp($line);
      my @line = split (/\s+/, $line);
      next unless ($#line);
      next if ($line =~ /^\s+/);
      unless ($line[0] =~ /_/) {
         $line[0] .= "_" . &maponegene($line[0]);
         $line = join ("  ", @line); # this just corrects $line if we change it
      }
      push (@dists, $line);
      my ($gene, $genome) = split (/_/, $line[0]);
      unless ($gene && $genome) {&niceexit("Can't parse $line in $file\n")}
      
      push (@genes, $gene);
      push (@genomes, $genome);
      $checkdupgenomes{$genome}++;
      }


   # now we loop through all the lines, and split them on white space. 
   # then we add each value to the pre-existing value in the matrix
   # note that because the genomes are represented as numbers we can just
   # use these numbers for the position in the matrix.
   # we are going to also count the number of times that we save each data
   # point for the final average.
   # Finally, we only do this in one direction because the input
   # matrices are complete (and identical) on both halves.
   
   # note that column zero of the matrix is empty (there is no genome 0)
   foreach my $z (0 .. $#dists) {
      my @line = split (/\s+/, $dists[$z]);
      unless ($#line == $#genomes) {
         my $x; foreach my $y (0 .. $#dists) {if ($dists[$y] eq $dists[$z]) {$x = $y}}
         &niceexit("PROBLEM WITH \n@line AND \n@genomes\n\nIN FILE: $file\n\nBECAUSE $#line AND $#genomes AT LINE $x\n");
         }
      my ($gene, $linegenome) = split (/_/, $line[0]);
      unless ($gene && $linegenome) {&niceexit("CAN'T PARSE @line SECOND TIME AROUND\n")}
      $linegenomecount{$linegenome}++;
      my @seengenome;
      foreach my $x (1 .. $#genomes) { #do this for all the genomes.
         next if ($x <= $z+1);
         # If we are padding the table with 100s where there is no match, we
         # need to convert the -1's to 100. Otherwise we will ignore it.
	 
	 
         if ($line[$x] == -1) {if ($pad) {$line[$x] = $penalty} else {next}}
	 next if ($line[$x] > $skip);
	 
  	 if ($args =~ /-l/) {
if ($args =~ /-c/) {print STDERR "For $x, $line[$x] has protein lengths $proteinlength{$genes[$x]} and $proteinlength{$gene} and becomes "}
           $line[$x] = $line[$x] * ($proteinlength{$genes[$x]}+$proteinlength{$gene});
	   unless ($args =~ /-a/) {$line[$x] = $line[$x]/2}
	   
if ($args =~ /-c/) {print STDERR " $line[$x]\n"}
	   if ($line[$x] > $maxmatch) {$maxmatch=$line[$x]}
	   if ($line[$x] < $minmatch) {$minmatch=$line[$x]}
	}
         #if it is itself, we want to make it zero. Otherwise, we'll save the protein numbers that match
         if ($genomes[$x] == $linegenome) {$line[$x] = '0.000'}
         else {
	    my $genematch;
            # save the protein matches, but I only want to save them one way around
            # to make it easier
            if ($gene <$genes[$x]) {$genematch = $gene.",".$genes[$x].";".$line[$x]}
               else {$genematch = $genes[$x].",".$gene.";".$line[$x]}
            # protein match is a two dimensional array where each element is an array.
            # but it is called with an array! 4 dimensions?
            ${$proteinmatches[$linegenome][$genomes[$x]]}{$genematch} =1;
            # gene matches is all the genes from $linegenome that match genome. This will
            # be used to calculate the penalty for ORFs that are missed.
            ${$genematches[$linegenome][$genomes[$x]]}{$gene} =1;
            }
         
	 # add the length if we need to.
#######         if ($args =~ /-l/ && $line[$x] > 0) 
	 if ($args =~ /-l/) {
	   #now save the data because the count is really the length not the number
	   $matrix[$linegenome][$genomes[$x]] += $line[$x];
           $newmatch[$linegenome][$genomes[$x]] += $line[$x];
	   
	   {
	   my $countadj;
	   if ($args =~ /-a/) {$countadj = ($proteinlength{$genes[$x]}+$proteinlength{$gene})}
	   else {$countadj = ($proteinlength{$genes[$x]}+$proteinlength{$gene})/2}
	   $count[$linegenome][$genomes[$x]] += $countadj;
	   $newmatchcount[$linegenome][$genomes[$x]] += $countadj;
	   }
	   
	   
         }
	 else {
           $matrix[$linegenome][$genomes[$x]] += $line[$x];
           $count[$linegenome][$genomes[$x]] ++;
         }
	 $seengenome[$linegenome][$genomes[$x]] ++;
            
         }

      # now we need to pad out all the missing genomes with 100's
      if ($pad) {
         foreach my $x (1 .. $nogenomes) {
           next if ($checkdupgenomes{$x});
           next if ($seengenome[$linegenome][$x]);
           if ($args =~ /-lp/) {
              $matrix[$linegenome][$x] += $penalty*$proteinlength{$gene};
              $count[$linegenome][$x] += $proteinlength{$gene};
	      $newmismatch[$linegenome][$x] += $penalty*$proteinlength{$gene};
	      $newmismatchcount[$linegenome][$x] += $proteinlength{$gene};
	      
           }
           else {
              $matrix[$linegenome][$x] += $penalty;
              $count[$linegenome][$x] ++;
           }
         }
      }
   }
}

print STDERR "\tDone\nSorting and calculating\n";
print STDERR "Minimum match was $minmatch and maximum match was $maxmatch\n";

my $genomeproteins;
{

# now we need to penalize genomes that have only a few macthes.
# we will go through gene matches for each pair in the matrix, and
# add a penalty based on the number of missing orfs.


	if ($pad) {

		if ($args =~ /-lp/) {$genomeproteins = &getallprots()}
		else {
			open (MISS, ">missing.seqs.txt") || &niceexit("Can't open missing.seqs.txt\n");
			print MISS "Original\t#ORFs\tCompared to\t# similar\t# different\n";
		}
		foreach my $y (0 .. $#genematches) {
			next unless (exists $noorfs{$y}); # this just checks we have orfs for genome $y
				foreach my $x (1 .. $#{$matrix[$y]}) {
					next unless (exists $noorfs{$x});
					next if ($y == $x);
					my @similar = keys %{$proteinmatches[$y][$x]};
					if ($y + $x ==19) {print "xxsimilar $y, $x -> ", $#similar+1, ":\n|", join ("|\n|", @similar), "|\n"}
# need to add a loop to get all proteins per genome, and then remove the ones we've seen
					if ($args =~ /-lp/) {
						my %found;
						foreach my $similar (@similar) {
							my ($genes, $trash) = split /;/, $similar;
							my ($gene1, $gene2) = split /,/, $genes;
							$found{$gene1}=$found{$gene2}=1;
						}

						foreach my $missedprot (@{${$genomeproteins}{$y}}) {
							next if (exists $found{$missedprot});
							$matrix[$y][$x] += $proteinlength{$missedprot}*$penalty;
							$count[$y][$x] += $proteinlength{$missedprot};


							$newmismatch[$y][$x] += $proteinlength{$missedprot}*$penalty;
							$newmismatchcount[$y][$x] += $proteinlength{$missedprot};
						}
					}
					else {
						my $difference = $noorfs{$y} - ($#similar+1);
						print MISS "$y\t$noorfs{$y}\t$x\t",$#similar+1, "\t$difference\n";
						next unless ($difference);
						$matrix[$y][$x] += ($penalty * $difference);
						$count[$y][$x] += $difference;
					}
				}
		}
	}
}

my %difference; my %genomedifference;

{

	my %seen;
# now we will average the matrix based on the count.
	foreach my $y (0 .. $#matrix) {
		next unless ($matrix[$y]);
		foreach my $x (1 .. $#{$matrix[$y]}) {
			next unless ($count[$y][$x] && $matrix[$y][$x]);
			my $temp = $x."+".$y; my $temp1 = $y."+".$x;
			next if ($seen{$temp} || $seen{$temp1});
			$seen{$temp} = $seen{$temp1} =1;

# because we are only looking at one half of the matrix (see above)
# we need to be sure that both halves are the same.
# this loop will take care of that.

			$matrix[$y][$x] = $matrix[$x][$y] = $matrix[$y][$x] + $matrix[$x][$y];
			$count[$y][$x] = $count[$x][$y] = $count[$y][$x] + $count[$x][$y];

			$matrix[$x][$y] = $matrix[$y][$x] = $matrix[$y][$x]/$count[$y][$x];


		}
	}

}
{
# we are going to output the matrix twice. This first loop will output the matrix with
# the replicates number, and the second loop will output the matrix alone with no replicates
# number. This is to test whether FITCH is breaking on the number of replicates.
	my $minmatch=100; my $maxmatch=1;
# now we have all the data, lets just print out the matrix
	open (OUT, ">$matrixfilename");
	print OUT $#matrix, "\n";
#foreach my $y (1 .. $#matrix) {print STDERR "\t$y"}
#print STDERR "\n";
	foreach my $y (1 .. $#matrix) {
		my $tempstring = "gnm".$y;
		if (length($tempstring) > 10) {print STDERR "$tempstring is too long\n"}
		my $spacestoadd = " " x (10 - length($tempstring));
		print OUT $tempstring,$spacestoadd;
		foreach my $x (1 .. $#matrix) {
			if ($y == $x)  {
				if ($args=~ /-l/) {
					my $total;
					foreach my $protein (@{${$genomeproteins}{$y}}) {$total+=$proteinlength{$protein}}
					print OUT "0 $total  ";
				}
				else {print OUT "0 $noorfs{$x}  "}
				next;
			}
			unless (defined $matrix[$y][$x]) {print OUT "$penalty 0  "; next}
			unless ($matrix[$y][$x]) {
				print OUT "0 ";
				if ($count[$y][$x]) {print OUT int($count[$y][$x]),"  "}
				else {print OUT "0  "}
				next;
			}
			if ($matrix[$y][$x] > $maxmatch) {$maxmatch=$matrix[$y][$x]}
			if ($matrix[$y][$x] < $minmatch) {$minmatch=$matrix[$y][$x]}
			print OUT $matrix[$y][$x], " ", int($count[$y][$x]), "  ";
		}
		print OUT "\n";
	}
	print STDERR "MATRIX: Minimum = $minmatch and MAXIMUM = $maxmatch\n";
	close OUT;
}

{
# output the matrix again, this time do not put out the replicates number
	open (OUT, ">$matrixfilename.nosubreplicates");
	print OUT $#matrix, "\n";
#foreach my $y (1 .. $#matrix) {print STDERR "\t$y"}
#print STDERR "\n";
	foreach my $y (1 .. $#matrix) {
		my $tempstring = "gnm".$y;
		if (length($tempstring) > 10) {print STDERR "$tempstring is too long\n"}
		my $spacestoadd = " " x (10 - length($tempstring));
		print OUT $tempstring,$spacestoadd;
		foreach my $x (1 .. $#matrix) {
			if ($y == $x)  {print OUT "0 "; next}
			unless (defined $matrix[$y][$x]) {print OUT "$penalty "; next}
			unless ($matrix[$y][$x]) {print OUT "0 "; next}
			print OUT $matrix[$y][$x], " ";
		}
		print OUT "\n";
	}
	close OUT;
}


if ($args =~ /-c/) {
	foreach my $key (sort {$difference{$a} <=> $difference{$b}} keys %difference) {print "$key difference: $difference{$key}\n"}
	foreach my $key (sort {$genomedifference{$b} <=> $genomedifference{$a}} keys %genomedifference) {print "genome$key difference: $genomedifference{$key}\n"}
}




if ($args=~ /-m/) {
   open (PROT, ">$dir.protein.matches") || &niceexit("Can't open $dir.protein.matches for writing\n");
   
   #print out all the protein matches
   foreach my $y (1 .. $nogenomes) {
   my $tempstring = "gnm".$y;
   if (length($tempstring) > 10) {print STDERR "$tempstring is too long\n"}
   my $spacestoadd = " " x (10 - length($tempstring));
   print PROT $tempstring,$spacestoadd, "\t";
   foreach my $x (1 .. $nogenomes) {
      unless (defined $proteinmatches[$y][$x]) {print PROT "\t"; next}
      unless ($proteinmatches[$y][$x]) {print PROT "\t"; next}
      my @allmatches = (keys %{$proteinmatches[$y][$x]}, keys %{$proteinmatches[$x][$y]});
      my %allmatches;
      @allmatches{@allmatches}=1;
      @allmatches = sort keys %allmatches;
      print PROT join (" ", sort @allmatches), "\t";
      }
   print PROT "\n";
   }
}




&niceexit(0);


sub readfile {
	open(IN, $genomedatafile) || die "can't open genomedatafile";
	while (<IN>)
	{
		next if (/^\#/);
		chomp;
		my @a=split /\t/;
		push @{$genomedata->{$a[0]}}, [$a[1], $a[2]];
	}
	close IN;
}

sub maponegene {
	&readfile() unless (defined $genomedata);
	my $gene=shift;
	foreach my $g (keys %$genomedata) {
		foreach my $tuple (@{$genomedata->{$g}}) {
			return $g if ($tuple->[0] eq $gene);
		}
	}
	return "";
}



sub getnoorfs {
	&readfile() unless (defined $genomedata);
	foreach my $g (keys %$genomedata) {
		$noorfs{$g}=scalar(@{$genomedata->{$g}});
	}
}

sub getprotlengths {
	&readfile() unless (defined $genomedata);
	local $| =1;
	my %length; my $total; my $count;
	print STDERR "Getting protein lengths ";

	foreach my $g (keys %$genomedata) {                
		foreach my $tuple (@{$genomedata->{$g}}) {
			$length{$tuple->[0]}=$tuple->[1];
			$total += $tuple->[1];
			$count++;
		}
	}

	print STDERR "Done\n";
	print STDERR "Total length found is $total for $count proteins, average is ", $total/$count, "\n";
	return \%length;
}



sub getallprots {
	my %genomeproteins;
	&readfile() unless (defined $genomedata);

	foreach my $g (keys %$genomedata) {
		foreach my $tuple (@{$genomedata->{$g}}) {
			push @{$genomeproteins{$g}}, $tuple->[0];
		}
	}
	return \%genomeproteins;
}

sub niceexit {
	my $reason = shift;
	if ($reason) {print STDERR $reason; exit(-1)}
	else {exit(0)}
}
