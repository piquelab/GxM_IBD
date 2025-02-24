library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)
library(annotables)
library(RColorBrewer)
library(Hmisc)

FDR <- 0.1

#overlap of ileum qnorm eqtls and INTACT results
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
#131





res_top <- res2 %>% filter(pid %in% intact_sig$Gene)

dim(res_top)
#51


write.table(res_top, file = "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD-eQTL_ileum_20PCs_eGenes_INTACT_IBD_overlap.txt", quote=F, sep="\t")

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
	png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/INTACT_eQTL_boxplots/IBD-eQTL__ileum_INTACT_IBD_overlap_significant_eGene_SNP_boxplot_20PCs_", gene_symbol,".png"))
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



