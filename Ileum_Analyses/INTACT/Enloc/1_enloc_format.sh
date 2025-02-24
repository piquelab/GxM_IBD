#!/bin/bash

cd $PWD

###
tissue=IBD_ileum
sbatch -q primary --mem=80G --time=1-23:00:00 -n 1 -N 1-1 --job-name=prepare_enloc_${tissue} --output=slurm_enloc.out --wrap "
   perl summarize_dap2enloc.pl -dir ./dap_rst_dir -vcf ../${tissue}.vcf.gz |gzip - > ${tissue}_fastenloc.eqtl.vcf.gz "
