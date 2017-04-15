#!/bin/bash
touch outfile; # this make sure that the outfile exists, so protdist asks for the correct name
/bin/echo -e "aligned/clustalw.$SGE_TASK_ID\nF\nprotdist/protdist.$SGE_TASK_ID\ny" | protdist

