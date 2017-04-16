#/bin/bash
perl ~/PhageProteomicTree/combineprotdists.pl protdist.fixed matrix $NUMGENOMES -d protein_genome_lengths.txt  -p 10 -l -lp > out 2> err
