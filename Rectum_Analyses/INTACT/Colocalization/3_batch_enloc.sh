#!/bin/bash

cd $PWD

### loop for traits
cat traits_ls.txt | \
while read trait;
do 
   ii=IBD_Rectum
   echo ${trait} ${ii}
   nsnp=`zcat ./gwas_PIP2/${trait}.pip.gz |wc -l`
   sbatch -q primary --mem=20G --time=2-10:00:00 -n 1 -N 1-1 --job-name=enloc_${trait}_${ii} --output=slurm_enloc_${trait}_${ii}.out --wrap "
   /wsu/home/ha/ha21/ha2164/Bin/fastenloc.static -eqtl ${ii}_fastenloc.eqtl.vcf.gz -go ./gwas_PIP2/${trait}.pip.gz -total_variants ${nsnp}  -prefix ./enloc_output/${trait}_${ii}"
   ###
done
