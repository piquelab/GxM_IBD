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

bcftools merge -l rnaseq_vcf.txt -Oz -o IBD_rnaseq_merged.vcf.gz


##renamed sample names
bcftools reheader -s IBD_eQTL_RNASeq_filtered_sample_names.txt -o merged_IBD_eQTL_rnaseq.vcf.gz IBD_rnaseq_merged.vcf.gz


bcftools view -i 'INFO/DP>100 & MAF>0.05' merged_IBD_eQTL_rnaseq.vcf.gz -Oz -o merged_IBD_eQTL_rnaseq.filter.vcf.gz



#intersect RNA and DNA genotypes and use those to filter RNA genotypes

module load samtools/1.11

bcftools query -f '%CHROM\t%POS0\t%END\t%ID\n' -o IBD_eQTL_RNA_merged.bed /wsu/home/groups/piquelab/IBD_eQTL/Genotype_QC/rnaseq_pileups/merged_IBD_eQTL_RNA_genotypes_filter.vcf.gz

module load bedtools

bedtools intersect -a IBD_eQTL_wxs_merged.bed -b IBD_eQTL_RNA_merged.bed > IBD_eQTL_common_vars.bed


module load samtools/1.11
bcftools view -R IBD_eQTL_common_vars.bed -Oz -o merged_IBD_eQTL_RNA_filtered_by_DNA_vars.vcf.gz /wsu/home/groups/piquelab/IBD_eQTL/Genotype_QC/rnaseq_pileups/merged_IBD_eQTL_RNA_genotypes_filter.vcf.gz


#merge DNA and RNA genotypes
module load samtools

bcftools merge /wsu/home/groups/piquelab/IBD_eQTL/Genotype_QC/merged_IBD_eQTL_RNA_filtered_by_DNA_vars.vcf.gz /wsu/home/groups/piquelab/IBD_eQTL/wxs_genotypes/IBD_eQTL_wxs_merged.vcf.gz -Oz -o IBD_eQTL_RNA_DNA_merged.vcf.gz

#Filtered RNA by DNA kinship
module load misc

plink2 --vcf merged_IBD_eQTL_RNA_filtered_by_DNA_vars.vcf.gz --maf 0.2 --make-king square  --make-king-table -out IBD_eQTL_filtered_RNA_kinship
#Merged RNA and DNA vcfs

plink2 --vcf IBD_eQTL_RNA_DNA_merged.vcf.gz --maf 0.20  --make-king square  --make-king-table --out IBD_eQTL_DNA_RNA_kinship
