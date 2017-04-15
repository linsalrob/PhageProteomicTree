# Phage Proteomic Tree for the Spounavirus analyses

The file spounavirus.tsv is the spreadsheet with all the spounavirus IDs in. Note that I added an extra column, the first column, with a unique number, that I will use while processing the sequences.


We prefered the RefSeq ID if possible, but otherwise we used the INSDC ID. To make a file with just the IDs we used this one line perl code:


```
perl -ne '@a=split /\t/; map {$a[$_] =~ s/^\s+//; $a[$_] =~ s/\s+$//} (2,3); 
          if ($a[2] eq "-") {print "$a[2]\n"} else {print "$a[3]\n"}' 
	     spounavirus.tsv > spounavirus_ids.txt
```

Next, we get all the genomes in GenBank format. For that, we use the command line [edirect utils from NCBI](https://www.ncbi.nlm.nih.gov/books/NBK179288/). (Note these weren't around in 1999!)

```
for GENOME in spounavirus_ids.txt;
   do echo $GENOME;
      esearch -db nucleotide -query $GENOME | efetch -format gb > genomes/$GENOME.gbk;
   done
```

Now we need to extract all the proteins from all the genomes, and create a list of names that we will use on the final tree. This list of names also includes an automatic abbreviation of the name that is generated. These names are not tested to be unique, but hopefully will be!


Then we create a list of [protein IDs, genome IDs, protein lengths] that we use when we combine everything later.


```
perl ../genbank_code/genbank2table.pl -t -p proteins.faa -i genome_names.txt genomes/*
perl ../protein_genome_length.pl proteins.faa  > protein_genome_length.txt

```

to test if the names are unique you can run this command and see if anything at the bottom of the list occurs more than once!

```
cut -f 3 -t$'\t' genome_names.txt | sort | uniq -c | sort -n 
```




