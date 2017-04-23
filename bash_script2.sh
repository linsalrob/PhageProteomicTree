#!/bin/bash

# when this is done concatenate everything
cat phage_proteins_blast/*blastp > phage_proteins.blastp

mkdir fastafiles
java -jar ~/PhageProteomicTree/ppt.jar phage_proteins.blastp 0.1 proteins.faa fastafiles

perl ~/PhageProteomicTree/make_id_maps.pl id.map protein_genome_length.txt genome_names.txt

export NUMFILES=$(ls fastafiles/ | sed -e 's/fasta.//' | sort -nr | head -n 1); echo "There are $NUMFILES files to process"
for DIR in sge_output aligned dnd protdist; do if [ ! -e $DIR ]; then mkdir $DIR; fi; done

