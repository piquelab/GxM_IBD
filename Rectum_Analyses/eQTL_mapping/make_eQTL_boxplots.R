library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)
library(annotables)
library(RColorBrewer)



##make rectum eQTL boxplots
# 1. txt file with the eQTL results to find out match SNPs with genes:
res <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/results/FastQTL_results_best_19GEPCs.txt",header=T, sep="\t", stringsAsFactors=FALSE)
res <- res[res$bqval<FDR,]
res <- res[order(res$bqval),]

# 2. residuals :
# 1. "bed" file with gene expression:
GE <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/IBD-eQTL_rectum_residuals_qnorm.txt", header=T, sep="\t", stringsAsFactors=FALSE)

cvf <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_eQTL_mapping_updated_genotypePCs_2_9_23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cvf <- as.data.frame(cvf)

cvf <- cvf %>% filter(SampleID %in% colnames(GE))
GE <- GE[,cvf$SampleID]

identical(cvf$SampleID, colnames(GE))

colnames(GE) <- cvf$SUBJECT_ID
# 3. txt file with genotypes:
gtbed <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/IBD-eQTL_PC19_signif_eGene_genotypes.txt",header=T, sep="\t", stringsAsFactors=FALSE)

gtbed$ref_length <- str_count(gtbed$ref)
gtbed$alt_length <- str_count(gtbed$alt)

gtbed <- gtbed[gtbed$ref_length == 1 & gtbed$alt_length == 1, ]

gtbed <- gtbed %>% select(-c(ref_length, alt_length))
gtbed <- gtbed[!duplicated(gtbed[,c("ID")]),]



##filter eQTL results for single variants
gtbed_snp_info <- gtbed %>% select(chr,pos,ID,ref,alt)
gtbed <- gtbed %>% select(-c(chr,pos,ID,ref,alt))



GE <- GE[,match(colnames(gtbed), colnames(GE))]


gtbed <- cbind(gtbed_snp_info, gtbed)

res_top <- head(res, n = 30)

for(i in 1:nrow(res_top)) {
	gene_id <- res_top$pid[i]
	gene_symbol <- grch38$symbol[grch38$ensgene == gene_id]
	snp_id <- res_top$sid[i]

	gt_data <- gtbed[gtbed$ID == snp_id,]
	ge_data <- GE[gene_id,]

	genotype_num <- t(gt_data)
	expression_data <- t(ge_data)

	genotype_num <- genotype_num[6:87,]

	df <- as.data.frame(cbind(genotype_num, expression_data))
	colnames(df) <- c("genotype_num", "expression_data")
	df$expression_data <- as.numeric(df$expression_data)
	df <- df %>% mutate(genotype = ifelse(genotype_num == "0/0", paste0(gt_data$ref, "/", gt_data$ref), 
											ifelse(genotype_num == "0/1", paste0(gt_data$ref, "/", gt_data$alt), 
												ifelse(genotype_num == "1/1", paste0(gt_data$alt, "/", gt_data$alt), "NA"))))
	mycolors <- c("#c25757", "#e9c61d","#3A68AE")
		names(mycolors) <- c(paste0(gt_data$ref, "/", gt_data$ref), paste0(gt_data$ref, "/", gt_data$alt),paste0(gt_data$alt, "/", gt_data$alt))

	df$genotype <- factor(df$genotype, levels=(c(paste0(gt_data$ref, "/", gt_data$ref), paste0(gt_data$ref, "/", gt_data$alt),paste0(gt_data$alt, "/", gt_data$alt))))

	png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/analysis/plots/IBD-eQTL_rectum_eQTL_boxplot_", gene_symbol, ".png"))
	p <- ggplot(df, aes(x=genotype, y=expression_data, fill = genotype)) + 
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width=0.25) + 
        scale_fill_manual(values=mycolors) +
        theme_classic() +
        ylab("Gene Expression Residuals") +
        ggtitle(paste(gene_symbol,snp_id)) +
        theme(legend.position="none", plot.title = element_text(hjust=0.5, size = rel(1.3)))
	print(p)
	dev.off()

}

