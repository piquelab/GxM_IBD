#!/bin/bash
#SBATCH -q primary
#SBATCH --mem=100G
#SBATCH --time=4-01:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=go7535@wayne.edu
#SBATCH -o output_%j.out
#SBATCH -e errors_%j.err

module load samtools

dbSnpFile=/nfs/rprdata/RefSnps/dbSNP144.hg38/snp144Common.bed.gz
genome=/wsu/home/groups/piquelab/pbmc_handls/ref_files/hg38.fa
bamFolder=/wsu/home/groups/piquelab/IBD_eQTL/bams

samtools mpileup -f ${genome} -l ${dbSnpFile} ${bamFolder}/${sample}_clean.bam -d 100000 -v -t DP,AD,ADF,ADR \
                | bcftools call - --ploidy GRCh38 -m \
                | bcftools view - -i 'FORMAT/DP>0' -Oz > ${sample}.vcf.gz
bcftools index ${sample}.vcf.gz 