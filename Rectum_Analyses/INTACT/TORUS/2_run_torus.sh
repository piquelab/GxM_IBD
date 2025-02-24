#!/bin/bash


###
### create output directory

outdir=./torus_output
if [ ! -d ${outdir} ]; then
   ##
   mkdir -p ${outdir}
fi


###
### submit jobs

sbatch -q primary --mem=80G --time=14-5:00:00 -n 1 -N 1-1 --job-name=torus_IBD_Rectum --output=slurm_torus_IBD_Rectum.output --wrap "
  module load misc;
  torus -d ./torus_input/IBD_Rectum.eQTL.txt.gz \
    -smap ./torus_input/zzz_snp.map.gz \
    -gmap ./torus_input/zzz_gene.map.gz \
    -est > ${outdir}/IBD_Rectum.est \
    -dump_prior ${outdir}/IBD_Rectum.dump.prior"

