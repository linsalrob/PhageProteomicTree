# Phage Proteomic Tree for the Spounavirus analyses

## Get the genomes and extract proteins

In this step, we will download all the genomes from genbank (note, they are in the [genomes](genomes/) directory already, so you can skip this step!), and extract protein sequences into a fasta file. We also get the information we need about the genomes.

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


## Pre-cluster the proteins based on blastp

In this step, we will use blastp to precluster the proteins. We also rename the proteins with a unique integer (as there an infinite number of integers!). In addition, some of the downstream steps have limited protein name lengths, and so using integers avoids any renaming issues.


Now we are going to identify our groups of proteins. We start by converting proteins.faa to a blast database. The makeblastdb code is part of [blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). The split_blast_queries_edwards_blastplus.pl is code that runs on our cluster using Sun Grid Engine to submit the blast to the cluster. You probably don't need to use the cluster, in which case you can just run a regular blast+ search to end up with the phage_proteins.blastp file.



```
makeblastdb -in proteins.faa -dbtype prot
perl ../split_blast_queries_edwards_blastplus.pl -f proteins.faa -n 100 -p blastp -db proteins.faa -d phage_proteins_blast -e 0.1 
# when this is done concatenate everything
cat phage_proteins_blast/*blastp > phage_proteins.blastp

```

We cluster all the proteins into groups based on their individual blast results. This is a greedy clustering, as in the next step we'll do an alignment and refine the results. I run this command on a machine with plenty of RAM!

```
mkdir fastafiles
java -jar ~/PhageProteomicTree/ppt.jar phage_proteins.blastp 0.1 proteins.faa fastafiles
```

This program writes out the fasta files in a directory (fastafiles). The fasta files have the clusters of proteins, but the proteins have been renamed with an integer for downstream analysis. The code also makes a file called `id.map` that has the original protein ids and the new protein ids (which are just integers).

We need to associate proteins with genomes, and genomes with IDs, so we make two map files, called genome_id.map and protein_genome_lengths.txt that you need for later steps.

```
perl ~/PhageProteomicTree/make_id_maps.pl id.map protein_genome_length.txt genome_names.txt
```

## Run clustalw and protdist

In this step, we will run clustalw and protdist. We use the cluster for this, but you probably don't have to since it is largely IO-bound anyway (it doesn't matter if you don't know what that means). If you need to download them, here are the links for [clustalw](http://www.clustal.org/download/current/) and [phylip](http://evolution.genetics.washington.edu/phylip/getme-new1.html) which contains protdist.

We figure out how many files we have, and we set that as an environment variable. We also set up some directories.

```
NUMFILES=$(ls fastafiles/ | sed -e 's/fasta.//' | sort -nr | head -n 1); echo "There are $NUMFILES files to process"
for DIR in sge_output aligned dnd protdist; do if [ ! -e $DIR ]; then mkdir $DIR; fi; done
```

Now comes the heavy lifting. We are going to run clustalw on all those sequence files and generate phylip output files. Then we're going to run protdist on all those output files to get our matrices. We use SGE for this, and submit to the cluster. The option -S says use bash as the shell, and the option -V says use your current environment variables, so that if clustalw or protdist are in your PATH they will be available on the nodes.

```
qsub -cwd  -S /bin/bash -V -o sge_output/ -e sge_output/ -t 1-$NUMFILES:1 ~/PhageProteomicTree/clustalw.sh
```

This gives us a job id, which we use as JID in this submission which holds until the alignments are done:

```
qsub -cwd -S /bin/bash -V -o sge_output/ -e sge_output/ -t 1-$NUMFILES:1 -hold_jid <previous job number> ~/PhageProteomicTree/protdist.sh
```























