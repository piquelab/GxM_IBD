#!/bin/bash/

condition=IBD_ileum

cat geneList_files.txt | \
while read geneFile;
do
echo ${geneFile}
sbatch --export=condition=${condition},geneFile=${geneFile} --output=slurm_assemble_${geneFile}.out --job-name=assemble_${condition}_${geneFile} 2_assemble_one.sh
sleep 0.5;
done


# chmod a+x ${condition}.assemble.cmd

