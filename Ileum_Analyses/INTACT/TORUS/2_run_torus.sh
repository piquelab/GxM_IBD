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
tissue=IBD_ileum
sbatch -q primary --mem=80G --time=14-5:00:00 -n 1 -N 1-1 --job-name=torus_${tissue} --output=slurm_torus_${tissue}.output --wrap "
  module load misc2;
  torus -d ./torus_input/${tissue}.eQTL.txt.gz \
    -smap ./torus_input/zzz_snp.map.gz \
    -gmap ./torus_input/zzz_gene.map.gz \
    -est > ${outdir}/${tissue}.est \
    -dump_prior ${outdir}/${tissue}.dump.prior"

