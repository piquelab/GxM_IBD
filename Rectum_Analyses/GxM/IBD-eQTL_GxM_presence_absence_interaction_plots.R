library(tidyverse)
library(qvalue)
library(ggplot2)
library(patchwork)
library(data.table)
library(stringr)
library(annotables)
library(Hmisc)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)

FDR <- 0.1




# 1. txt file with the significant interaction eQTLs
signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/signif_interactions_list_uncorrected.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("trait", "k,", "Interactions")

signif <- signif %>% filter(Interactions > 0)
dim(signif)
#14
sum(signif$Interactions)
#19

traits <- signif$trait

# 2. residuals :
# 1. "bed" file with gene expression:
GE <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/IBD-eQTL_rectum_residuals_qnorm.txt", header=T, sep="\t", stringsAsFactors=FALSE)

cvf <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_eQTL_mapping_updated_genotypePCs_2_9_23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cvf <- as.data.frame(cvf)

cvf <- cvf %>% filter(SampleID %in% colnames(GE))
GE <- GE[,cvf$SampleID]

identical(cvf$SampleID, colnames(GE))

colnames(GE) <- cvf$SUBJECT_ID

wxs_cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD_eQTL_WXS_metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
wxs_cov <- as.data.frame(wxs_cov)

wxs_rectum <- wxs_cov %>% filter(submitted_subject_id %in% cvf$SUBJECT_ID) %>% arrange(Run)

GE <- GE[,wxs_rectum$submitted_subject_id]
# 3. txt file with genotypes:
gtbed <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/IBD-eQTL_PC19_signif_eGene_genotypes.txt",header=T, sep="\t", stringsAsFactors=FALSE)

gtbed$ref_length <- str_count(gtbed$ref)
gtbed$alt_length <- str_count(gtbed$alt)

gtbed <- gtbed[gtbed$ref_length == 1 & gtbed$alt_length == 1, ]

gtbed <- gtbed %>% select(-c(ref_length, alt_length))

#dosages
# 2. txt file with dosages:
dosbed <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/IBD-eQTL_rectum_PC19_signif_eGene_dosages.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
dosbed <- dosbed[!duplicated(dosbed[,c("ID")]),]
rownames(dosbed) <- dosbed[,3]
dos <- dosbed[,-c(1:3)]

# 4. microbial abundance values:
# 4. microbial abundance values:
cv <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/dichotomous_microbes/IBD-eQTL_Rectum_covariates_dichotomous_microbe_4-24-23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cv <- as.data.frame(cv)

#subset dosages and GE to only samples we have microbial data for

cv <-cv[match(colnames(dos), cv$SUBJECT_ID),]
rownames(cv) <- cv$SUBJECT_ID

#subset dosages and GE to only samples we have microbial data for
gtbed_snp_info <- gtbed[,1:5]
gtbed <- gtbed[,6:91]

gtbed <- gtbed[,colnames(gtbed) %in% cv$SUBJECT_ID]

GE <- GE[,colnames(GE) %in% cv$SUBJECT_ID]
dos <- dos[,colnames(dos) %in% cv$SUBJECT_ID]

cv <-cv[match(colnames(gtbed), cv$SUBJECT_ID),]
rownames(cv) <- cv$SUBJECT_ID


GE <- GE[,match(colnames(gtbed), colnames(GE))]
dos <- dos[,match(colnames(gtbed), colnames(dos))]

gtbed <- cbind(gtbed_snp_info, gtbed)

##make boxplots for signficant trait-gene-snp interactions

for(t in 1:length(traits)) {
	trait <- traits[t]
	results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxE_presence_absence_19PCs_", trait, ".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

	cv_trait <- cv[,c(trait)]

	
	#filter for significant interactions using  FDR
	res_sig <- results %>% filter(interaction_qval<FDR)
	for(i in 1:nrow(res_sig)) {
		gene_id <- res_sig$pid[i]
		gene_symbol <- grch38$symbol[grch38$ensgene == gene_id]
		snp_id <- res_sig$sid[i]

		gt_data <- gtbed[gtbed$ID == snp_id,]
		ge_data <- GE[gene_id,]
		dos_data <- dos[snp_id, ]

		genotype_num <- t(gt_data)
		expression_data <- t(ge_data)
		dosages <- t(dos_data)

		genotype_num <- genotype_num[6:75,]

		df <- as.data.frame(cbind(genotype_num, expression_data, dosages))
		colnames(df) <- c("genotype_num", "expression_data", "dosage")
		df$expression_data <- as.numeric(df$expression_data)
		df <- df %>% mutate(genotype = ifelse(genotype_num == "0/0", paste0(gt_data$ref, "/", gt_data$ref), 
												ifelse(genotype_num == "0/1", paste0(gt_data$ref, "/", gt_data$alt), 
													ifelse(genotype_num == "1/1", paste0(gt_data$alt, "/", gt_data$alt), "NA"))))
		df$var <- cv_trait

		df <- df %>% mutate(dichotomous = ifelse(var == 0, "Absent", "Present"))

		model <- lm(expression_data~dosages*cv_trait)
		Intercept <- summary(model)$coefficients[1,1]
		dosage_beta <- summary(model)$coefficients[2,1]
		metagene_beta <- summary(model)$coefficients[3,1]
		interaction_beta  <- summary(model)$coefficients[4,1]
		in0 <- Intercept
		in1 <- Intercept+dosage_beta
		in2 <- Intercept+2*dosage_beta
		slop0 <- metagene_beta
		slop1 <- metagene_beta+interaction_beta
		slop2 <- metagene_beta+2*interaction_beta

		mycolors <- c("#c25757", "#e9c61d","#3A68AE")
		names(mycolors) <- c(paste0(gt_data$ref, "/", gt_data$ref), paste0(gt_data$ref, "/", gt_data$alt),paste0(gt_data$alt, "/", gt_data$alt))

		df$genotype <- factor(df$genotype, levels=(c(paste0(gt_data$ref, "/", gt_data$ref), paste0(gt_data$ref, "/", gt_data$alt),paste0(gt_data$alt, "/", gt_data$alt))))
		

		png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/plots/interaction_plots/Interaction_eQTL_plots_", trait,"_",gene_symbol, "_GxE_presence_absence_19PCs.png"))
		p <- ggplot(df, aes(x=genotype, y= expression_data, fill = genotype)) +
			geom_boxplot() +
			geom_smooth(method='lm', se = FALSE, color = "black", aes(group=1)) +
			geom_point(color="black", size=1.0, alpha=1.0) +
			theme_classic() +
			scale_fill_manual(values=mycolors) +
			xlab(trait) +
			ylab(paste(gene_symbol, " Gene Expression Residuals")) +
			ggtitle(paste(gene_symbol, snp_id)) +
			theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), axis.title.x = element_text(size = 10)) +
			facet_wrap(~dichotomous)
		print(p)
		dev.off()
	} 
}



#make heatmap of interaction effects

#create results file with all GxE results for dichotomous test
vars <- read.table("/wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_protein_coding_residuals/GxE_dichotomous_19PCs/1-per-var-list.txt", sep="\t", comment="=", blank.lines.skip=T)

colnames(vars) <- "Var"

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxE_presence_absence_19PCs_", vars$Var[1], ".txt")
res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
res$var <- vars$Var[1]


res_all <- data.frame(matrix(ncol = ncol(res), nrow = 0))
colnames(res_all) <- colnames(res)

for(i in 1:length(vars$Var)) {
	var <- vars$Var[i]
	fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxE_presence_absence_19PCs_", var, ".txt")
	res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
	res$var <- var

	anno <- grch38 %>% select(ensgene,symbol) %>% filter(grch38$biotype == "protein_coding") %>% unique() %>% as.data.frame()
	anno <- anno %>% filter(ensgene %in% res$pid) %>% dplyr::rename(pid = ensgene)
	res <- merge(res, anno, by="pid")

	res_all <- rbind(res_all, res)

}

write.table(res_all, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxE_presence_absence_all_results_19PCs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

all_results_presence_absence <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxE_presence_absence_all_results_19PCs.txt", sep="\t", stringsAsFactors=F, header=T, comment="")

res_sig <- all_results_presence_absence %>% filter(interaction_qval < FDR)

gene_taxa_sig <- res_sig %>% filter(var %in% traits) %>% select(var, symbol) %>% arrange(var)

write.table(gene_taxa_sig, file="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxM_presence_absence_19PCs_sig_genes_and_taxa_summary.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

genes <- res_sig %>% select(pid, symbol) %>% unique()
dim(genes)
#11 2

write.table(genes, file="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxM_presence_absence_19PCs_sig_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

traits
# [1] "Bacteria.Firmicutes.Negativicutes.Selenomonadales.Veillonellaceae.Veillonella"
 [2] "Bacteria.Proteobacteria.Betaproteobacteria.Burkholderiales.Alcaligenaceae.Sutterella"
 [3] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae..Eubacterium..hallii.group"
 [4] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Subdoligranulum"
 [5] "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Porphyromonadaceae.Parabacteroides"
 [6] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae.ND3007.group"
 [7] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Fusicatenibacter"
 [8] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
 [9] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Erysipelotrichaceae.UCG.003"
[10] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae..Eubacterium..eligens.group"


microbe_short <- c("Veillonella", "Sutterella", "Eubacterium.hallii.group", "Subdoligranulum", "Parabacteroides", "Lachnospiraceae.ND3007.group", 
	"Fusicatenibacter", "Lachnospira", "Erysipelotrichaceae.UCG.003", "Eubacterium.eligens.group")

gene_list <- c()

res_list <- list()
IE_cols <- c()
normalized_IE_cols <- c()

for(i in 1:length(traits)) {
	trait <- traits[i]
	fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxE_presence_absence_19PCs_", trait, ".txt")
	res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
	res$microbe <- trait
	res$microbe_short <- microbe_short[i]
	res$ie_zscore <- res$interaction_effect / res$interaction_SE
	res <- res %>% dplyr::rename(!!paste0("taxa",i,"_GxM_IE") := interaction_effect)
	res <- res %>% dplyr::rename(!!paste0("taxa",i,"_GxM_IE_zscore") := ie_zscore)
	res_list[[i]] <- res
	IE_cols <- c(IE_cols, paste0("taxa",i,"_GxM_IE"))
	normalized_IE_cols <- c(normalized_IE_cols, paste0("taxa",i,"_GxM_IE_zscore"))

	res_sig <- res %>% filter(interaction_qval < FDR)
	gene_list <- c(gene_list, res_sig$pid)

}

taxa_id_table <- cbind(traits, IE_cols, normalized_IE_cols)
write.table(taxa_id_table, file="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/GxE_presence_absence_17PCs/plots/interaction_plots/significant_taxa_id_table.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

for(i in 1:length(traits)) {
	res_list[[i]] <- res_list[[i]] %>% filter(pid %in% gene_list)
}

res_list <- res_list %>% reduce(full_join, by="pid")

res_sc <- res_list[,c(IE_cols)]

colnames(res_sc) <- microbe_short

res_sc_cor <- rcorr(as.matrix(res_sc), type = "spearman")

res_sc_cor_val <- res_sc_cor$r
res_sc_cor_p <- res_sc_cor$P

#check they have same dimensions
dim(res_sc_cor_val) == dim(res_sc_cor_p)

##make the nonsignificant correlations NAs
nonsig <- which(res_sc_cor_p > 0.05, arr.ind = T)

res_cor_val_sig <- res_sc_cor_val

for(i in 1:nrow(nonsig)) {
    row_num <- nonsig[i,1]
    col_num <- nonsig[i,2]
    res_cor_val_sig[row_num, col_num] <- 0
}
#make heatmap of correlations, make 0s gray - i.e. the non signficiant correlations

library(pheatmap)
pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/plots/interaction_plots/GxM_dichotomous_interaction_effect_correlation_heatmap.pdf") 
p <- pheatmap(res_sc_cor_val, fontsize = 9)
print(p) 
dev.off()

pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/plots/interaction_plots/GxM_dichotomous_interaction_effect_signficant_correlation_heatmap.pdf") 
p <- pheatmap(res_cor_val_sig, fontsize = 9)
print(p) 
dev.off()

res_zscore <- res_list[,c(normalized_IE_cols)]

colnames(res_zscore) <- microbe_short

res_zscore_cor <- rcorr(as.matrix(res_zscore), type = "spearman")

res_zscore_cor_val <- res_zscore_cor$r
res_zscore_cor_p <- res_zscore_cor$P

#check they have same dimensions
dim(res_zscore_cor_val) == dim(res_zscore_cor_p)

##make the nonsignificant correlations NAs
nonsig <- which(res_zscore_cor_p > 0.05, arr.ind = T)

res_zscore_cor_val_sig <- res_zscore_cor_val

for(i in 1:nrow(nonsig)) {
    row_num <- nonsig[i,1]
    col_num <- nonsig[i,2]
    res_zscore_cor_val_sig[row_num, col_num] <- 0
}
#make heatmap of correlations, make 0s gray - i.e. the non signficiant correlations

bk2 <- c(0.01,seq(0.01,1,by=0.1))
my_palette <- c("gray90", "gray90",
              c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/plots/interaction_plots/GxM_dichotomous_normalized_IE_correlation_heatmap_significant.pdf") 
p <- pheatmap(res_zscore_cor_val_sig,color = my_palette, fontsize = 9)
print(p) 
dev.off()

