#!/bin/bash

rm -f `find protdist -size 0`
perl ~/PhageProteomicTree/rewrite_protdists.pl protdist protdist.fixed
export NUMGENOMES=$(wc -l genome_names.txt | sed -e 's/\s\+genome_names.txt//')
qsub -cwd -S /bin/bash -V -v NUMGENOMES=$NUMGENOMES -o sge_output -e sge_output ~/PhageProteomicTree/matrix.sh
