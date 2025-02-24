#!/bin/bash/


cd $PWD
 
cat geneList_files.txt |
while read geneFile; do
   ##
   sbatch -q express -p erprp --mem=4G --time=1-00:00:00 -N 1-1 -n 1 --job-name=topPIP_${geneFile} --output=slurm_topPIP_${geneFile}.output --wrap "
   module load R;
   R CMD BATCH --no-save --no-retore '--args ${geneFile}' 5_topPIP_one.R topPIP_${geneFile}.Rout"
   echo ${geneFile}
   sleep 0.5;
done 


###
# cd 5_summary.outs
# cat ./IBD_Rectum/dap*.txt > IBD_Rectum_PIP.txt
# gzip IBD_Rectum_PIP.txt
