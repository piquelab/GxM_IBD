#!/bin/bash/



cat geneList_files.txt | \
while read geneFile; do
##
echo ${geneFile} 
sbatch --export=geneFile=${geneFile} --job-name=parse_${geneFile} --output=slurm_parse_${geneFile}.out 4_parse_one.sh
##sleep 1;

done





