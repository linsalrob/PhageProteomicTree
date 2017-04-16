#/bin/bash
echo "There are $NUMGENOMES genomes to process"
perl ~/PhageProteomicTree/combineprotdists.pl protdist.fixed matrix $NUMGENOMES -d protein_genome_lengths.txt  -p 10 -l -lp 
