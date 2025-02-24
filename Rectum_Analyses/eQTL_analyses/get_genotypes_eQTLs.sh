module load samtools/1.11

bcftools view -i'ID=@/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/analysis/PC19_significant_topeeQTL_snps.txt' /wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_bi-allelic_SNPs/IBD_eQTL_rectum_DNA_genotypes_filtered_SNPs.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' > PC19_eGene_genotypes.txt