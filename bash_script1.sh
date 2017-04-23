#!/bin/bash

perl ~/PhageProteomicTree/genbank_code/genbank2table.pl -t -p proteins.faa -i genome_names.txt genomes/*
perl ~/PhageProteomicTree/protein_genome_length.pl proteins.faa  > protein_genome_length.txt

cut -f 3 -t$'\t' genome_names.txt | sort | uniq -c | sort -n 

makeblastdb -in proteins.faa -dbtype prot
perl ~/PhageProteomicTree/split_blast_queries_edwards_blastplus.pl -f proteins.faa -n 100 -p blastp -db proteins.faa -d phage_proteins_blast -e 0.1 

