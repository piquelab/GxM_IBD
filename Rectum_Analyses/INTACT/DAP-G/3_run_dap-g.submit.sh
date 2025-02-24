#!/bin/bash/

condition=IBD_Rectum
cat geneList_files.txt | \
while read geneFile;
do
echo ${condition} ${geneFile}
sbatch --export=condition=${condition},geneFile=${geneFile} --output=slurm_dap_${geneFile}.out --job-name=dap_${condition}_${geneFile} 3_run_dap-g.one.sh
sleep 1;
done
