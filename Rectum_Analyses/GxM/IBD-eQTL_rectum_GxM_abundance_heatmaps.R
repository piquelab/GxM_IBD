library(tidyverse)
library(ggplot2)
library(igraph)
library(data.table)
library(stringr)
library(annotables)
library(RColorBrewer)
library(ggraph)
library(Hmisc)
library(irlba)
library(pheatmap)

FDR <- 0.1

#heatmap
vars <- read.table("/wsu/home/groups/piquelab/IBD_eQTL/covariates/abundance_model_all_microbes_1-per-var-list.txt", sep="\t", comment="=", blank.lines.skip=T)

colnames(vars) <- "Var"

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_", vars$Var[1], ".txt")
res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
res$var <- vars$Var[1]


res_all <- data.frame(matrix(ncol = ncol(res), nrow = 0))
colnames(res_all) <- colnames(res)


for(i in 1:length(vars$Var)) {
  var <- vars$Var[i]
  fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_", var, ".txt")
  res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
  res$var <- var
  res_all <- rbind(res_all, res)
}

anno <- grch38 %>% filter(ensgene %in% unique(res_all$pid)) %>% select(ensgene, symbol) %>% dplyr::rename(pid = ensgene)

res_all <- merge(res_all, anno, by = "pid")

write.table(res_all, file="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_all_results.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

redundant_microbes <- c("Bacteria.Firmicutes.Clostridia.Clostridiales", "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales", "Bacteria.Proteobacteria.Deltaproteobacteria", "Bacteria.Proteobacteria.Deltaproteobacteria.Desulfovibrionales.Desulfovibrionaceae")

res_all <- res_all %>% filter(!var %in% redundant_microbes)

res_sig <- res_all %>% filter(interaction_qval < FDR)



gene_taxa_sig <- res_sig %>% select(var, symbol) %>% arrange(var)
dim(gene_taxa_sig)

 # 25  2


write.table(gene_taxa_sig, file="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxM_abundance_19PCs_sig_genes_and_taxa_summary.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

genes <- res_sig %>% select(pid, symbol) %>% unique()
dim(genes)
#23 2

write.table(genes, file="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxM_abundance_19PCs_sig_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

gene_list <- genes$pid

traits <- unique(gene_taxa_sig$var)


traits

#  [1] "Bacteria.Firmicutes.Clostridia.Clostridiales.Family.XIII.Family.XIII.AD3011.group"
#  [2] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae"           
#  [3] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
#  [4] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceaencultured"  
#  [5] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae"           
#  [6] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Butyricicoccus"
#  [7] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Ruminiclostridium.6"
#  [8] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceaencultured"  
#  [9] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae"
# [10] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Holdemanella"
# [11] "Bacteria.Firmicutes.Negativicutes.Selenomonadales.Acidaminococcaceae.Phascolarctobacterium"
# [12] "Bacteria.Proteobacteria.Deltaproteobacteria.Desulfovibrionales.Desulfovibrionaceae.Bilophila"
# [13] "Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae"
# [14] "Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae.Haemophilus"


microbe_short <- c("Clostridiales.Family.XIII.AD3011", "Lachnospiraceae", "Lachnospira", "Lachnospiraceaencultured", "Ruminococcaceae",  "Butyricicoccus","Ruminiclostridium.6", 
  "Ruminococcaceaencultured", "Erysipelotrichaceae","Holdemanella", "Phascolarctobacterium", "Bilophila", "Pasteurellaceae", "Haemophilus")

res_list <- list()
IE_cols <- c()
normalized_IE_cols <- c()

for(i in 1:length(traits)) {
	trait <- traits[i]
	fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_", trait,".txt")
	res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
	res <- res %>% filter(pid %in% gene_list)
	res$microbe <- trait
  res$microbe_short <- microbe_short[i]
	res$ie_zscore <- res$interaction_effect / res$interaction_SE
	res <- res %>% dplyr::rename(!!paste0("taxa",i,"_GxM_IE") := interaction_effect)
	res <- res %>% dplyr::rename(!!paste0("taxa",i,"_GxM_IE_zscore") := ie_zscore)
	res_list[[i]] <- res
	IE_cols <- c(IE_cols, paste0("taxa",i,"_GxM_IE"))
	normalized_IE_cols <- c(normalized_IE_cols, paste0("taxa",i,"_GxM_IE_zscore"))

}

taxa_id_table <- cbind(traits, IE_cols, normalized_IE_cols)
write.table(taxa_id_table, file="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/significant_taxa_id_table.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

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
bk1 <- c(seq(-1, -0.01, by = 0.1), -0.001)  # Negative values (down to slightly below 0)
bk2 <- c(0.001, seq(0.01, 1, by = 0.1))    # Positive values (from slightly above 0)
bk <- c(bk1, bk2)  # Combine break limits

# Create custom palette: blue for negative, gray for zero, red for positive
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1) - 1),
                "gray90", "gray90",  # Assign gray for near-zero correlations
                colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2) - 1))



library(pheatmap)
pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_interaction_effect_heatmap.pdf") 
p <- pheatmap(res_sc_cor_val, color = my_palette,  # Use custom palette
              breaks = bk,         # Define unique custom breaks
              fontsize = 9)
print(p) 
dev.off()

pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_interaction_effect_signficant_heatmap.pdf") 
p <- pheatmap(res_cor_val_sig, color = my_palette,  # Use custom palette
              breaks = bk,         # Define unique custom breaks
              fontsize = 9)
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

#make heatmap of correlations, make 0s gray - i.e. the non signficiant correlations
bk1 <- c(seq(-1, -0.01, by = 0.1), -0.001)  # Negative values (down to slightly below 0)
bk2 <- c(0.001, seq(0.01, 1, by = 0.1))    # Positive values (from slightly above 0)
bk <- c(bk1, bk2)  # Combine break limits

# Create custom palette: blue for negative, gray for zero, red for positive
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1) - 1),
                "gray90", "gray90",  # Assign gray for near-zero correlations
                colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2) - 1))



pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_normalized_IE_heatmap.pdf") 
p <- pheatmap(res_zscore_cor_val, color = my_palette,  # Use custom palette
              breaks = bk,         # Define unique custom breaks
              fontsize = 9)
print(p) 
dev.off()

pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_normalized_IE_heatmap_significant.pdf") 
p <- pheatmap(res_zscore_cor_val_sig, color = my_palette,  # Use custom palette
              breaks = bk,         # Define unique custom breaks
              fontsize = 9)
print(p) 
dev.off()

dendrogram <- p$tree_row

# Cut the dendrogram into two clusters
clusters <- cutree(dendrogram, k = 2)

# Count the number of observations in each cluster
cluster_counts <- table(clusters)

# Sort the cluster counts in descending order
cluster_counts <- sort(cluster_counts, decreasing = TRUE)

# Extract the two major clusters
major_clusters <- names(cluster_counts)[1:3]

cluster_df <- data.frame(names = names(clusters), cluster = clusters, row.names = NULL)


 PCs <- prcomp_irlba(t(res_zscore), n=5)


 summary(PCs)
 mypcs <- as.data.frame(PCs$x)
 rownames(mypcs) <- cluster_df$names


PCvar <- data.frame(summary(PCs)$importance)[2,]
PCvar <- t(PCvar)

meta <- cbind(cluster_df, mypcs)

meta <- meta %>% mutate(cluster_id = ifelse(cluster == 1, "1", ifelse(cluster == 2, "2", "3")))


percentVar <- round(PCvar * 100)

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/IBD-eQTL_rectum_GxM_abundance_normalized_interaction_effects_PCS1and2_Plot.pdf")
pdf(fname)
pca <- ggplot(meta, aes(PC1, PC2, color=cluster_id)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("GxM Abundance normalized IE PCs") +
  theme_classic() +
  coord_fixed()
print(pca) 
dev.off()

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/IBD-eQTL_rectum_GxM_abundance_normalized_interaction_effects_taxa_labeled_PCS1and2_Plot.pdf")
pdf(fname)
pca <- ggplot(meta, aes(PC1, PC2, color=names, shape=cluster_id)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  coord_fixed()
print(pca) 
dev.off()

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/IBD-eQTL_rectum_GxM_abundance_normalized_interaction_effects_taxa_labeled_no_legend_PCS1and2_Plot.pdf")
pdf(fname)
pca <- ggplot(meta, aes(PC1, PC2, color=names, shape=cluster_id)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none")
print(pca) 
dev.off()




cv <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD-eQTL_Rectum_microbiome_all_taxa_covariates_6-7-23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cv <- as.data.frame(cv)


#subset columns for microbe abundance

traits <- traits[!traits %in% redundant_microbes]

cv_microbe <- cv[,c(traits)]

colnames(cv_microbe)



colnames(cv_microbe) <- microbe_short
res_sc_cor <- rcorr(as.matrix(cv_microbe), type = "spearman")

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
bk1 <- c(seq(-0.7,-0.1,by=0.1),-0.001)
bk2 <- c(0.001,seq(0.1,1,by=0.1))
bk <- c(bk1,bk2)  #combine the break cs for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
              "gray90", "gray90",
              c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))


png("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/GxM_microbe_abundance_correlation_heatmap.png") 
p <- pheatmap(res_sc_cor_val, color = my_palette,fontsize = 9)
print(p) 
dev.off()

png("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/GxM_microbe_abundance_signficant_correlation_heatmap.png") 
p <- pheatmap(res_cor_val_sig, color = my_palette, fontsize = 9)
print(p) 
dev.off()


microbe_PCs <- prcomp_irlba(t(cv_microbe), n=5)


 summary(microbe_PCs)
 mypcs <- as.data.frame(microbe_PCs$x)
 rownames(mypcs) <- colnames(cv_microbe)
 mypcs$microbe_short <- rownames(mypcs)


PCvar <- data.frame(summary(microbe_PCs)$importance)[2,]
PCvar <- t(PCvar)


meta <- cbind(cluster_df, mypcs)

meta <- meta %>% mutate(cluster_id = ifelse(cluster == 1, "1", ifelse(cluster == 2, "2", "3")))


percentVar <- round(PCvar * 100)

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/IBD-eQTL_rectum_GxM_CLR_abundance_PCS1and2_Plot.pdf")
pdf(fname)
pca <- ggplot(meta, aes(PC1, PC2, color=cluster_id)) +
  geom_point(size=5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  ggtitle("GxM taxa CLR abundance PCs") +
  theme_classic() +
  coord_fixed()
print(pca) 
dev.off()

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/IBD-eQTL_rectum_GxM_CLR_abundance_taxa_labeled_PCS1and2_Plot.pdf")
pdf(fname)
pca <- ggplot(meta, aes(PC1, PC2, color=names, shape=cluster_id)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  coord_fixed()
print(pca) 
dev.off()

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/IBD-eQTL_rectum_GxM_CLR_abundance_taxa_labeled_no_legend_PCS1and2_Plot.pdf")
pdf(fname)
pca <- ggplot(meta, aes(PC1, PC2, color=names, shape=cluster_id)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none")
print(pca) 
dev.off()


#make interaction plots for examples of genes in taxa in opposite cluster

goi <- "ABO"

#import results of all taxa
res_abundance <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_all_results.txt", sep="\t", stringsAsFactors=F, header=T, comment="")


#filter for taxa in cluster 2 and gene of interest - check to see if interaction p-value is significant

cluster_df$trait <- traits

cluster2_microbes <- cluster_df %>% filter(cluster == 2) %>% select(trait) %>% unlist()

res_filtered <- res_abundance[res_abundance$var %in% cluster2_microbes & res_abundance$symbol == "ABO", ]

res_filtered_sig <- res_filtered %>% filter(interaction_pval < 0.05)

res_filtered_sig 
#  pid         sid Intercept_pval  dosage_pval metagene_pval
# 1 ENSG00000175164 rs950529388   2.074844e-06 5.064696e-08    0.06743361
# 2 ENSG00000175164 rs950529388   2.449323e-06 5.801868e-08    0.06505786
#   interaction_pval Intercept_effect dosage_effect metagene_effect
# 1       0.04512704        0.9699322    -0.7721474       0.1452771
# 2       0.04439638        0.9585101    -0.7659044       0.1483472
#   interaction_effect Intercept_SE dosage_SE metagene_SE interaction_SE
# 1         -0.1108242    0.1863891 0.1255284  0.07813263     0.05426610
# 2         -0.1116129    0.1857493 0.1252044  0.07907011     0.05446001
#   Intercept_qval  dosage_qval metagene_qval interaction_qval
# 1   4.263511e-05 3.758190e-07     0.8467551        0.6359065
# 2   4.929018e-05 4.514283e-07     0.8325439        0.6138858
#                                                                                      var
# 1 Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae.Haemophilus
# 2             Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae
#   symbol
# 1    ABO
# 2    ABO


#make interaction plots


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
cv <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD-eQTL_Rectum_microbiome_all_taxa_covariates_6-7-23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
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

for(t in 1:nrow(res_filtered_sig)) {
  trait <- res_filtered_sig$var[t]
  results <- read.table(paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_", trait,".txt"), sep="\t", header=TRUE, stringsAsFactors=FALSE)

  cv_trait <- cv[,c(trait)]

    gene_id <- res_filtered_sig$pid[t]
    gene_symbol <- res_filtered_sig$symbol[t]
    snp_id <- res_filtered_sig$sid[t]

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
    png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/Example_for_manuscript_Abundance_GxM_19PCs_Interaction_eQTL_Scatterplot_", trait,"_",gene_symbol, ".png"))
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
    png(paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/Example_for_manuscript_Abundance_GxM_19PCs_Interaction_eQTL_tertile_boxplot_", trait, "_",gene_symbol,".png"))
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


#print summary of lm models for gene taxa
library(tidyverse)
library(ggplot2)
library(igraph)
library(data.table)
library(stringr)

res_abundance <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_all_results.txt", sep="\t", stringsAsFactors=F, header=T, comment="")

traits_of_interest <- c("Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Butyricicoccus", "Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae.Haemophilus")

res_filtered <- res_abundance[res_abundance$var %in% traits_of_interest & res_abundance$symbol == "ABO", ]

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
cv <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD-eQTL_Rectum_microbiome_all_taxa_covariates_6-7-23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
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

for(t in 1:nrow(res_filtered)) {
  trait <- res_filtered$var[t]
cv_trait <- cv[,c(trait)]

    gene_id <- res_filtered$pid[t]
    gene_symbol <- res_filtered$symbol[t]
    snp_id <- res_filtered$sid[t]

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

    model <- lm(expression_data~dosages*cv_trait)
    print(trait)
    print(summary(model))
}


cv_trait1 <- cv[,c((traits_of_interest[1]))]
cv_trait2 <- cv[,c((traits_of_interest[2]))]


gene_id <- res_filtered$pid[1]
    gene_symbol <- res_filtered$symbol[1]
    snp_id <- res_filtered$sid[1]

    gt_data <- gtbed[gtbed$ID == snp_id,]
    ge_data <- GE[gene_id,]
    dos_data <- dos[snp_id, ]

    genotype_num <- t(gt_data)
    expression_data <- t(ge_data)
    dosages <- t(dos_data)

    genotype_num <- genotype_num[6:75,]

  

    model <- lm(expression_data~ dosages + cv_trait1 + cv_trait2 + dosages:cv_trait1 + dosages:cv_trait2)
       print(summary(model))