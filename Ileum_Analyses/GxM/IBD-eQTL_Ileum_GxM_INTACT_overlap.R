##Script: IBD-eQTL_Ileum_GxM_INTACT_overlap.R
##Author: Shreya Nirmalan
##Last Modified: 10-25-2024
##Purpose: to overlap GxM QTLs with eGenes thar are associated with IBD, CD and UC risk based on integration with INTACT results for Ileum eQTL data

library(tidyverse)
library(ggplot2)
library(igraph)
library(data.table)
library(stringr)
library(annotables)
library(RColorBrewer)
library(ggraph)
library(qvalue)

dataDir <- "/rs/rs_grp_ibdeqtl/IBD_shreya_2024/IBD_ileum_2024-10-22/4_TWAS_smr/INTACT_output/eqtl_topPIP"

traits <- c("IBD.EUR.Crohns_Disease", "IBD.EUR.Inflammatory_Bowel_Disease", "IBD.EUR.Ulcerative_Colitis")

resList <- list()
#read in files and extract signficiant risk associated eGenes

for(i in 1:length(traits)) {
	trait <- traits[i]
	fname <- paste0(dataDir, "/", trait, "_IBD_ileum_intact.txt")
	res <- read.table(fname, stringsAsFactors=F, header=T, comment="")

			res <- res %>% select(Gene, GRCP, GLCP, zscore, LFDR, FDR)
			res$Trait <- trait


			resList[[i]] <- res

		}


intact_all <- data.frame(matrix(ncol = ncol(resList[[1]]), nrow = 0))
        colnames(intact_all) <- colnames(resList[[1]])

        for(n in 1:length(resList)) {
            intact_all <- rbind(intact_all, resList[[n]])
        }


intact_sig <- intact_all %>% filter(FDR < 0.1)
length(intact_sig$Gene) 
#209
length(unique(intact_sig$Gene))

intact_sig <- intact_sig %>% dplyr::rename(ensgene = Gene)

anno <- grch38 %>% select(ensgene,symbol) %>% filter(grch38$biotype == "protein_coding") %>% unique() %>% as.data.frame()
anno <- anno %>% filter(ensgene %in% intact_sig$ensgene)
	intact_sig <- merge(intact_sig, anno, by="ensgene")

write.csv(intact_sig, "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD_eQTL_ileum_significant_INTACT_genes.csv")

#import abundance GxE results to make overlap
all_results_abundance <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_all_results.txt", sep="\t", stringsAsFactors=F, header=T, comment="")
FDR <- 0.1
res_sig <- all_results_abundance %>% filter(interaction_qval < FDR)

res_sig <- res_sig %>% dplyr::rename(microbe = var)


#make microbe_gene column
res_sig$microbe_gene <- paste0(res_sig$microbe,"_",res_sig$symbol)


#check overlap with INTACT
table(res_sig$pid %in% intact_sig$ensgene)


# FALSE  TRUE
#    23     1
#check overlap with all intact results
table(res_sig$pid %in% intact_all$Gene)

# TRUE
#   24

#merge GxM with INTACT

intact_all_gene <- intact_all %>% dplyr::rename(pid = Gene)

intact_gxm_abundance <- merge(intact_all_gene, res_sig, by="pid")

intact_gxm_abundance <- intact_gxm_abundance %>% select(pid, symbol, microbe, Trait, zscore,GLCP, FDR, interaction_effect, interaction_qval)

write.csv(intact_gxm_abundance, "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD_eQTL_intact_gxm_abundance.csv")


##check if INTACT genes have GxM

res_intact <- all_results_abundance %>% filter(pid %in% intact_sig$ensgene)

#filter p-val significance
res_intact_sig <- res_intact %>% filter(interaction_pval < 0.05)

dim(res_intact_sig)
#[1] ] 341  20

head(res_intact_sig)

#filter intact for microebs with GxM
res_intact_sig_gxm <- res_intact_sig %>% filter(var %in% res_sig$microbe)
dim(res_intact_sig_gxm)
#34 20

length(unique(res_intact_sig_gxm$pid))
#22 genes 
length(unique(res_intact_sig_gxm$var))
#12

res_intact_sig_gxm <- res_intact_sig_gxm %>% dplyr::rename(microbe = var)
intact_sig_gene <- intact_sig %>% dplyr::rename(pid = ensgene)
intact_sig_gxm_abundance <- merge(intact_sig_gene, res_intact_sig_gxm, by="pid")

intact_sig_gxm_abundance <- intact_sig_gxm_abundance %>% select(pid, symbol.x, microbe, Trait, zscore, GLCP, FDR, interaction_effect, interaction_pval, interaction_qval)

write.csv(intact_sig_gxm_abundance, "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD_eQTL_intact_sig_gxm_abundance_pval_overlap.csv")


##P/A GxM
res <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_presence_absence/GxE_results/Ileum_GxE_presence_absence_all_results_20PCs.txt", sep="\t", stringsAsFactors=F, header=T, comment="")

res_sig <- res %>% filter(interaction_qval < 0.1)

res_sig <- res_sig %>% dplyr::rename(microbe = var)



table(unique(res$pid) %in% unique(intact_sig$ensgene))

# FALSE  TRUE
#  3261    51


table(res_sig$pid %in% intact_sig$ensgene)
# FALSE
#    10
table(res_sig$pid %in% intact_all$Gene)
# TRUE    10

intact_gxm_dichotomous <- merge(intact_all_gene, res_sig, by="pid")

intact_gxm_dichotomous <- intact_gxm_dichotomous %>% select(pid, symbol, microbe, Trait, zscore, GLCP, FDR, interaction_effect, interaction_qval)

write.csv(intact_gxm_dichotomous, "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD_eQTL_intact_gxm_dichotomous.csv")


res_intact <- res %>% filter(pid %in% intact_sig$ensgene)

#filter p-val significance
res_intact_sig <- res_intact %>% filter(interaction_pval < 0.05)

dim(res_intact_sig)
# [1] 150  20


head(res_intact_sig)

#filter intact for microebs with GxM
res_intact_sig_gxm <- res_intact_sig %>% filter(var %in% res_sig$microbe)
dim(res_intact_sig_gxm)
#] 23 20


length(unique(res_intact_sig_gxm$pid))
#15 genes 
length(unique(res_intact_sig_gxm$var))
#8

res_intact_sig_gxm <- res_intact_sig_gxm %>% dplyr::rename(microbe = var)
intact_sig_gene <- intact_sig %>% dplyr::rename(pid = ensgene)
intact_sig_gxm_abundance <- merge(intact_sig_gene, res_intact_sig_gxm, by="pid")

intact_sig_gxm_abundance <- intact_sig_gxm_abundance %>% select(pid, symbol.x, microbe, Trait, zscore, GLCP, FDR, interaction_effect, interaction_pval, interaction_qval)

write.csv(intact_sig_gxm_abundance, "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD_eQTL_intact_sig_gxm_presence_absence_pval_overlap.csv")




##re-calculate GxM padj for INTACT genes for abundance 

intact_genes <- unique(intact_sig$ensgene)

microbes <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/1-per-var-list.txt", sep="\t", stringsAsFactors=F, comment="")

microbe_list <- microbes$V1

res_summary <- data.frame(matrix(ncol=2,nrow=0))

colnames(res_summary) <- c("microbe", "sigGxM")

for(i in 1:length(microbe_list)) {

	microbe <- microbe_list[i]

	fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/IBD-eQTL_Ileum_GxE_abundance_20PCs_", microbe, ".txt")
  results <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
	
	res_intact <- results %>% filter(pid %in% intact_genes)

	res_intact$new_interaction_qval <- p.adjust(res_intact$interaction_pval, method = "bonferroni")

	write.table(res_intact, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/INTACT_ileum_sig_genes_GxE_abundance_20PCs_", microbe, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

	numGxM <- nrow(res_intact %>% filter(new_interaction_qval < 0.1))

	vec <- c(microbe, numGxM)

	res_summary[nrow(res_summary) + 1, ] <- vec


}

write.table(res_summary, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/INTACT_ileum_sig_genes_GxE_abundance_20Cs_summary.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

signif <- res_summary %>% filter(sigGxM > 0)
dim(signif)
#9 taxa
signif <- signif %>% distinct(microbe, .keep_all = TRUE)



traits <- signif$microbe

#make GxM abundance plots

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
	results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/INTACT_ileum_sig_genes_GxE_abundance_20PCs_", trait, ".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

	cv_trait <- cv[,c(trait)]


	#filter for significant interactions using  FDR
	res_sig <- results %>% filter(new_interaction_qval<0.1)
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
		png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/plots/INTACT_Ileum_Abundance_GxM_20PCs_Interaction_eQTL_Scatterplot_", trait,"_",gene_symbol, ".png"))
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
	results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/INTACT_ileum_sig_genes_GxE_abundance_20PCs_", trait, ".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

	cv_trait <- cv[,c(trait)]

	
	#filter for significant interactions using  FDR
	res_sig <- results %>% filter(new_interaction_qval<0.1)
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
		png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/plots/INTACT_Ileum_Abundance_GxM_20PCs_Interaction_eQTL_tertile_boxplot_", trait, "_",gene_symbol,".png"))
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



#P/A gxm INTACT


microbes <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/presence-absence-1-per-var-list.txt", sep="\t", stringsAsFactors=F, comment="")

microbe_list <- microbes$V1

res_summary <- data.frame(matrix(ncol=2,nrow=0))

colnames(res_summary) <- c("microbe", "sigGxM")

for(i in 1:length(microbe_list)) {

	microbe <- microbe_list[i]

	fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_presence_absence/GxE_results/IBD-eQTL_Ileum_GxE_presence_absence_20PCs_", microbe, ".txt")
	results <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
	
	res_intact <- results %>% filter(pid %in% intact_genes)

	res_intact$new_interaction_qval <- p.adjust(res_intact$interaction_pval, method = "bonferroni")

	write.table(res_intact, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/presence_absence/INTACT_ileum_sig_genes_GxE_presence_absence_20PCs_", microbe, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

	numGxM <- nrow(res_intact %>% filter(new_interaction_qval < 0.1))

	vec <- c(microbe, numGxM)

	res_summary[nrow(res_summary) + 1, ] <- vec


}

write.table(res_summary, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/presence_absence/INTACT_ileum_sig_genes_GxE_presence_absence_20PCs_summary.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

signif <- res_summary %>% filter(sigGxM > 0)
dim(signif)
#10
signif <- signif %>% distinct(microbe, .keep_all = TRUE)


traits <- signif$microbe

FDR <- 0.1

# 1. "bed" file with gene expression:
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
cv <- read_delim("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/IBD-eQTL_Ileum_covariates_with_presence_absence_microbial_data_10-2-2024.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
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
	results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/presence_absence/INTACT_ileum_sig_genes_GxE_presence_absence_20PCs_", trait, ".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

	cv_trait <- cv[,c(trait)]

	#filter for significant interactions using  FDR
	res_sig <- results %>% filter(new_interaction_qval<FDR)
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
		png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/presence_absence/plots/INTACT_Ileum_Interaction_eQTL_plots_", trait,"_", gene_symbol, "_GxE_presence_absence_20PCs.png"))
		
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


#make summary table of INTACT genes wth GxM - the sign of interaction effect and direction of risk

#abundance GxM
signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/INTACT_ileum_sig_genes_GxE_abundance_20Cs_summary.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("microbe", "sigGxM")

signif <- signif[-1,]

signif <- signif %>% filter(sigGxM > 0)
dim(signif)
#12 taxa
signif <- signif %>% distinct(microbe, .keep_all = TRUE)


traits <- signif$microbe



results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/INTACT_ileum_sig_genes_GxE_abundance_20PCs_", traits[1], ".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
results$microbe <- traits[1]

INTACT_sig_abundance_GxM <- data.frame(matrix(ncol = ncol(results), nrow = 0))
colnames(INTACT_sig_abundance_GxM) <- colnames(results)


for(t in 1:length(traits)) {
	trait <- traits[t]
	results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/abundance/INTACT_ileum_sig_genes_GxE_abundance_20PCs_", trait, ".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

	res_sig <- results %>% filter(new_interaction_qval < 0.1)
	res_sig$microbe <- trait

	INTACT_sig_abundance_GxM <- rbind(INTACT_sig_abundance_GxM, res_sig)

}


INTACT_sig_abundance_GxM$GxM <- "Abundance GxM"

#PA GxM


signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/presence_absence/INTACT_ileum_sig_genes_GxE_presence_absence_20PCs_summary.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("microbe", "sigGxM")

signif <- signif[-1,]

signif <- signif %>% filter(sigGxM > 0)
dim(signif)
#12 taxa
signif <- signif %>% distinct(microbe, .keep_all = TRUE)


traits <- signif$microbe


results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/presence_absence/INTACT_ileum_sig_genes_GxE_presence_absence_20PCs_", traits[1], ".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)
results$microbe <- traits[1]

INTACT_sig_pa_GxM <- data.frame(matrix(ncol = ncol(results), nrow = 0))
colnames(INTACT_sig_pa_GxM) <- colnames(results)


for(t in 1:length(traits)) {
	trait <- traits[t]
	results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/GxM_INTACT_results/presence_absence/INTACT_ileum_sig_genes_GxE_presence_absence_20PCs_", trait, ".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

	res_sig <- results %>% filter(new_interaction_qval < 0.1)
	res_sig$microbe <- trait

	INTACT_sig_pa_GxM <- rbind(INTACT_sig_pa_GxM, res_sig)

}

INTACT_sig_pa_GxM$GxM <- "P/A GxM"

#get INTACT gene info

INTACT_sig_GxM <- rbind(INTACT_sig_abundance_GxM, INTACT_sig_pa_GxM)

#merge w INTACT results
INTACT_sig <- read.csv("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD_eQTL_ileum_significant_INTACT_genes.csv")

INTACT_sig <- INTACT_sig %>% select(ensgene, symbol, GLCP, zscore, FDR, Trait) %>% dplyr::rename(pid=ensgene)

INTACT_sig_GxM <- merge(INTACT_sig, INTACT_sig_GxM, by="pid")

write.table(INTACT_sig_GxM, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD-eQTL_Ileum_INTACT_sig_genes_with_GxM_full_results.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

INTACT_sig_GxM_short <- INTACT_sig_GxM %>% select(pid, symbol, zscore, FDR, Trait, microbe, interaction_effect, interaction_pval, new_interaction_qval)

write.table(INTACT_sig_GxM_short, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD-eQTL_Ileum_INTACT_sig_genes_with_GxM_condensed_results.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write.csv(INTACT_sig_GxM_short, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD-eQTL_Ileum_INTACT_sig_genes_with_GxM_condensed_results.csv"))


#make boxplots for INTACT genes with GxM

#get list of genes

INTACT_GxM <- read.csv("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD-eQTL_Ileum_INTACT_sig_genes_with_GxM_condensed_results.csv")

# 1. txt file with the eQTL results to find out match SNPs with genes:
res <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/permutations/results/FastQTL_results_best_20GEPCs.txt",header=T, sep="\t", stringsAsFactors=FALSE)
res <- res[res$bqval<FDR,]
res <- res[order(res$bqval),]

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



##filter eQTL results for single variants
res2 <- res %>% filter(sid %in% gtbed$ID)


GE_data <- GE[,match(colnames(gtbed)[6:87], colnames(GE))]
res_top <- res2 %>% filter(pid %in% INTACT_GxM$pid)

dim(res_top)


for(i in 1:nrow(res_top)) {
	gene_id <- res_top$pid[i]
	gene_symbol <- grch38$symbol[grch38$ensgene == gene_id][1]
	snp_id <- res_top$sid[i]

	gt_data <- gtbed[gtbed$ID == snp_id,]
	if(nrow(gt_data > 1)) {
		gt_data <- gt_data[1,]
	}
	ge_data <- GE_data[gene_id,]

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
	png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/INTACT_eQTL_boxplots/IBD-eQTL_ileum_INTACT_IBD_with GxM_genes_eQTL_boxplot_", gene_symbol,".png"))
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
