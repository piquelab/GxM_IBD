#!/bin/bash/

### create directory
outdir=./fastQTL.output/
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi

###
### fastqtl
tissue=IBD_Rectum
for i in {1..30}; do
     sbatch -q express -p erprp --mem=10G --time=24:00:00 -N 1-1 -n 1 --job-name=${tissue}_chunk${i} --out=slurm_fastQTL_chunk${i}.output --wrap "
     module load misc2;
     fastQTL --vcf ../${tissue}.vcf.gz --bed ../${tissue}.bed.gz --out ${outdir}/${tissue}_nominals_chunk${i}.txt.gz --window 1e6 --chunk ${i} 30 \
     --cov ../${tissue}_PC1-19.covariates.txt "
     sleep 0.5;
done
