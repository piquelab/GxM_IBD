#!/bin/bash

cd $PWD

###
sbatch -q primary --mem=80G --time=1-23:00:00 -n 1 -N 1-1 --job-name=prepare_enloc --output=slurm_enloc.out --wrap "
   perl summarize_dap2enloc.pl -dir ./dap_rst_dir -vcf ../IBD_Rectum.vcf.gz |gzip - > IBD_Rectum_fastenloc.eqtl.vcf.gz "
