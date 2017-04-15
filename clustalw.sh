#!/bin/bash
clustalw -infile=fastafiles/fasta.$SGE_TASK_ID -outfile=aligned/clustalw.$SGE_TASK_ID -output=phylip -newtree=dnd/tree.$SGE_TASK_ID.dnd -align
