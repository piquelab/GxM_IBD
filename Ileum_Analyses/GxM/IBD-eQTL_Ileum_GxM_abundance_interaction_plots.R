library(tidyverse)
library(ggplot2)
library(igraph)
library(data.table)
library(stringr)
library(annotables)
library(RColorBrewer)
library(ggraph)
library(Hmisc)


##figure 4 abundance GxM scatterplots. tertile boxplots and heatmap
FDR <- 0.1

# 1. txt file with the significant interaction eQTLs
signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/IBD-eQTL_Ileum_GxM_abundance_interaction_summary.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("trait", "Interactions")

signif <- signif %>% filter(Interactions > 0)
dim(signif)
#14 taxa
signif <- signif %>% distinct(trait, .keep_all = TRUE)
sum(signif$Interactions)
#24 interactions

table(signif$Interactions > 1)
#5 taxa with more than 1 interaction

traits <- signif$trait




# 2. residuals :
# 1. "bed" file with gene expression:
GE <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/IBD-eQTL_ileum_residuals_q_normed.txt", header=T, sep="\t", stringsAsFactors=FALSE)

cvf <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_eQTL_mapping_updated_genotypePCs_2_9_23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cvf <- as.data.frame(cvf)

cvf <- cvf %>% filter(SampleID %in% colnames(GE))
GE <- GE[,cvf$SampleID]

identical(cvf$SampleID, colnames(GE))

colnames(GE) <- cvf$SUBJECT_ID

wxs_cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD_eQTL_WXS_metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
wxs_cov <- as.data.frame(wxs_cov)

wxs_ileum <- wxs_cov %>% filter(submitted_subject_id %in% cvf$SUBJECT_ID) %>% arrange(Run)

GE <- GE[,wxs_ileum$submitted_subject_id]

# 3. txt file with genotypes:
gtbed <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/permutations/IBD-eQTL_Ileum_PC20_signif_eGene_genotypes.txt",header=T, sep="\t", stringsAsFactors=FALSE)

gtbed$ref_length <- str_count(gtbed$ref)
gtbed$alt_length <- str_count(gtbed$alt)

gtbed <- gtbed[gtbed$ref_length == 1 & gtbed$alt_length == 1, ]

gtbed <- gtbed %>% select(-c(ref_length, alt_length))

#dosages
# 2. txt file with dosages:
dosbed <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/permutations/IBD-eQTL_ileum_PC20_signif_dosages.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
dosbed <- dosbed[!duplicated(dosbed[,c("ID")]),]
rownames(dosbed) <- dosbed[,3]
dos <- dosbed[,-c(1:3)]

# 4. microbial abundance values:
cv <- read_delim("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/IBD-eQTL_Ileum_covariates_with_microbial_data_10-2-2024.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cv <- as.data.frame(cv)

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

##make scatterplots for signficant trait-gene-snp interactions

for(t in 1:length(traits)) {
	trait <- traits[t]
	results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/IBD-eQTL_Ileum_GxE_abundance_20PCs_", trait,".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

	cv_trait <- cv[,c(trait)]


	#filter for significant interactions using  FDRR
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

		genotype_num <- genotype_num[6:74,]

		df <- as.data.frame(cbind(genotype_num, expression_data, dosages))
		colnames(df) <- c("genotype_num", "expression_data", "dosage")
		df$expression_data <- as.numeric(df$expression_data)
		df <- df %>% mutate(genotype = ifelse(genotype_num == "0/0", paste0(gt_data$ref, "/", gt_data$ref), 
												ifelse(genotype_num == "0/1", paste0(gt_data$ref, "/", gt_data$alt), 
													ifelse(genotype_num == "1/1", paste0(gt_data$alt, "/", gt_data$alt), "NA"))))
		df$var <- cv_trait

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
		png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/scatterplots/Ileum_Abundance_GxM_20PCs_Interaction_eQTL_Scatterplot_", trait,"_",gene_symbol, ".png"))
		p <- ggplot(df, aes(x=var, y= expression_data, color = genotype)) +
			geom_point(alpha=0.7)+
			geom_abline(intercept=in0,slope=slop0,color="#c25757", size=1)+
			geom_abline(intercept=in1,slope=slop1,color="#e9c61d", size=1)+
			geom_abline(intercept=in2,slope=slop2,color="#3A68AE", size=1)+
			theme_classic() +
			scale_colour_manual(values=mycolors) +
			xlab(trait) +
			ylab(paste(gene_symbol, "Gene Expression Residuals")) +
			ggtitle(paste(snp_id)) +
			theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), axis.title.x = element_text(size = rel(0.7)))

		print(p)
		dev.off()
	}
	
}

for(t in 1:length(traits)) {
	trait <- traits[t]
	results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/IBD-eQTL_Ileum_GxE_abundance_20PCs_", trait,".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

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

		genotype_num <- genotype_num[6:74,]

		df <- as.data.frame(cbind(genotype_num, expression_data, dosages))
		colnames(df) <- c("genotype_num", "expression_data", "dosage")
		df$expression_data <- as.numeric(df$expression_data)
		df <- df %>% mutate(genotype = ifelse(genotype_num == "0/0", paste0(gt_data$ref, "/", gt_data$ref), 
												ifelse(genotype_num == "0/1", paste0(gt_data$ref, "/", gt_data$alt), 
													ifelse(genotype_num == "1/1", paste0(gt_data$alt, "/", gt_data$alt), "NA"))))
		df$var <- cv_trait

		t <- quantile(df$var, probs = seq(0,1, (1/3)))

		df <- df %>% mutate(abundance = ifelse(var < t[2], "Low Abundance",
								ifelse(var > t[3], "High Abundance", "Medium Abundance")))
		df$abundance <- factor(df$abundance, levels = c("Low Abundance", "Medium Abundance", "High Abundance"))

		mycolors <- c("#c25757", "#e9c61d","#3A68AE")
		names(mycolors) <- c(paste0(gt_data$ref, "/", gt_data$ref), paste0(gt_data$ref, "/", gt_data$alt),paste0(gt_data$alt, "/", gt_data$alt))

		df$genotype <- factor(df$genotype, levels=(c(paste0(gt_data$ref, "/", gt_data$ref), paste0(gt_data$ref, "/", gt_data$alt),paste0(gt_data$alt, "/", gt_data$alt))))
		png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/tertile_boxplots/Ileum_Abundance_GxM_20PCs_Interaction_eQTL_tertile_boxplot_", trait, "_",gene_symbol,".png"))
		p <- ggplot(df, aes(x=genotype, y= expression_data, fill = genotype)) +
			geom_boxplot() +
			geom_point(color="black", size=1.0, alpha=1.0) +
			theme_classic() +
			scale_fill_manual(values=mycolors) +
			xlab(trait) +
			ylab(paste(gene_symbol, " Gene Expression Residuals")) +
			ggtitle(paste(gene_symbol, snp_id)) +
			theme(plot.title = element_text(hjust=0.5, size = rel(1.3)), axis.title.x = element_text(size = 10)) +
			facet_wrap(~abundance)
		print(p)
		dev.off()
	}
	
}


