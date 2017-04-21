#!/bin/bash

mkdir neighbor
cp matrix.nosubreplicates neighbor/infile
cd neighbor
echo -e "j\n133\ny"  | neighbor
cp outtree ../raw.tree
cd ..
perl ~/PhageProteomicTree/rename_tree_leaves.pl genome_id.map raw.tree > renamed_full.tree

