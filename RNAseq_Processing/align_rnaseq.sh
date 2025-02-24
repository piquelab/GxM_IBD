#!/bin/bash
#SBATCH -q express -p erprp
#SBATCH --mem=100G
#SBATCH --time=4-01:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=go7535@wayne.edu
#SBATCH -o output_%j.out
#SBATCH -e errors_%j.err

module load samtools

genomeindex=/nfs/rprdata/Anthony/data/HISAT2_Index/grch38_snp_tran/genome_snp_tran
###Align Reads###
module load hisat2
hisat2 -p 8 -x ${genomeindex} -1 /nfs/rprdata/download/HMP2-IBD/rnaseq/fastq/${var}_1.fastq.gz \
                              -2 /nfs/rprdata/download/HMP2-IBD/rnaseq/fastq/${var}_2.fastq.gz \
      2> ${var}_aligned.bam.e | samtools view -b1 - > ${var}_aligned.bam



###Sort Reads###
samtools sort -@ 4 -T tmp_${var}_aligned.bam -o ${var}_sorted.bam ${var}_aligned.bam
samtools index ${var}_sorted.bam
samtools view -c ${var}_sorted.bam > ${var}_sorted_count.txt


###Quality Filter###
samtools view -b1 -q10 ${var}_sorted.bam > ${var}_quality.bam
samtools index ${var}_quality.bam
samtools view -c ${var}_quality.bam > ${var}_quality_count.txt


###Deduplication###
samtools rmdup ${var}_quality.bam ${var}_clean.bam
samtools index ${var}_clean.bam
samtools view -c ${var}_clean.bam > ${var}_clean_count.txt

