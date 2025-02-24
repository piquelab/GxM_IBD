#!/bin/bash/

condition=IBD_ileum
cat geneList_files.txt | \
while read geneFile;
do
echo ${condition} ${geneFile}
sbatch --export=condition=${condition},geneFile=${geneFile} --output=slurm_dap_${geneFile}.out --job-name=dap_${condition}_${geneFile} 3_run_dap-g.one.sh
sleep 0.5;
done
