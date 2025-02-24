for j in $(seq 1 30); do
    sbatch -q primary -N1-1 -n 2 --mem=12G -t 10000 --job-name=$i$j --wrap "module load misc; fastQTL --vcf /wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_bi-allelic_SNPs/IBD_eQTL_rectum_DNA_genotypes_filtered_SNPs.vcf.gz --bed /rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/IBD-eQTL_Rectum_covariate_corrected_voom-q-norm.bed.gz --permute 1000 10000 --window 1e5 --out output/PC1-0.permutations.chunk$j.txt.gz --chunk $j 30"
sleep 1
done