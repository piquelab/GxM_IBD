#!/bin/bash/


condition=IBD_ileum
echo ${condition} 
sbatch -q express -p erprp -N 1-1 -n 1 --mem=5G --time=1:00:00 --job-name=${condition} --wrap "
      perl batch_process.pl -e ../${condition}.bed.gz \
       -g  ../${condition}.vcf.gz \
       -c ../${condition}_PC1-20.covariates.txt -t ${condition}"
