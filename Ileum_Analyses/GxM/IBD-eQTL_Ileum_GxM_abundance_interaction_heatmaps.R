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
vars <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/1-per-var-list.txt", sep="\t", comment="=", blank.lines.skip=T)

colnames(vars) <- "Var"

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/IBD-eQTL_Ileum_GxE_abundance_20PCs_", vars$Var[1], ".txt")
res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
res$var <- vars$Var[1]


res_all <- data.frame(matrix(ncol = ncol(res), nrow = 0))
colnames(res_all) <- colnames(res)


for(i in 1:length(vars$Var)) {
  var <- vars$Var[i]
  fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/IBD-eQTL_Ileum_GxE_abundance_20PCs_", var, ".txt")
  res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
  res$var <- var
  res_all <- rbind(res_all, res)
}

anno <- grch38 %>% filter(ensgene %in% unique(res_all$pid)) %>% select(ensgene, symbol) %>% dplyr::rename(pid = ensgene)

res_all <- merge(res_all, anno, by = "pid")

write.table(res_all, file="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_all_results.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



res_sig <- res_all %>% filter(interaction_qval < FDR)



gene_taxa_sig <- res_sig %>% select(var, symbol) %>% arrange(var)
dim(gene_taxa_sig)

 # 24  2


write.table(gene_taxa_sig, file="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_sig_genes_and_taxa_summary.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

genes <- res_sig %>% select(pid, symbol) %>% unique()
dim(genes)
#23 2

write.table(genes, file="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_sig_genes.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

gene_list <- genes$pid

traits <- unique(gene_taxa_sig$var)


traits

# [1] "Bacteria.Actinobacteria.Actinobacteria.Micrococcales"                                  
#  [2] "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella_9"          
#  [3] "Bacteria.Firmicutes.Bacilli.Bacillales.Family_XI.Gemella"                              
#  [4] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Anaerostipes"             
#  [5] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Coprococcus_1"            
#  [6] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnoclostridium"        
#  [7] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_008"  
#  [8] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_010"  
#  [9] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceaencultured"                 
# [10] "Bacteria.Firmicutes.Clostridia.Clostridiales.Peptostreptococcaceae.Intestinibacter"    
# [11] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Clostridium_innocuum_group"
# [12] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Erysipelotrichaceae_UCG_003"
# [13] "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae"      
# [14] "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae.Escherichia_Shigella"

redundant_microbes <- c("Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae")

traits <- traits[!traits %in% redundant_microbes]

microbe_short <- c("Micrococcales", "Prevotella_9", "Bacillales.Family_XI.Gemella", "Anaerostipes", "Coprococcus_1", "Lachnoclostridium", "Lachnospiraceae_UCG_008", "Lachnospiraceae_UCG_010", "Lachnospiraceaencultured", "Intestinibacter", "Clostridium_innocuum_group", "Erysipelotrichaceae_UCG_003","Escherichia_Shigella")

res_list <- list()
IE_cols <- c()
normalized_IE_cols <- c()

for(i in 1:length(traits)) {
	trait <- traits[i]
	fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/IBD-eQTL_Ileum_GxE_abundance_20PCs_", trait,".txt")
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
write.table(taxa_id_table, file="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/significant_taxa_id_table.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

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
pdf("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_interaction_effect_heatmap.pdf") 
p <- pheatmap(res_sc_cor_val, color = my_palette, breaks = bk, fontsize = 9)
print(p) 
dev.off()

pdf("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_interaction_effect_signficant_heatmap.pdf") 
p <- pheatmap(res_cor_val_sig, color = my_palette, breaks = bk, fontsize = 9)
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
bk1 <- c(seq(-1, -0.01, by = 0.1), -0.001)  # Negative values (down to slightly below 0)
bk2 <- c(0.001, seq(0.01, 1, by = 0.1))    # Positive values (from slightly above 0)
bk <- c(bk1, bk2)  # Combine break limits

# Create custom palette: blue for negative, gray for zero, red for positive
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1) - 1),
                "gray90", "gray90",  # Assign gray for near-zero correlations
                colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2) - 1))



pdf("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_normalized_IE_heatmap.pdf") 
p <- pheatmap(res_zscore_cor_val, color = my_palette, breaks = bk, fontsize = 9)
print(p) 
dev.off()

pdf("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_normalized_IE_heatmap_significant.pdf") 
p <- pheatmap(res_zscore_cor_val_sig, color = my_palette, breaks = bk, fontsize = 9)
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

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/IBD-eQTL_ileum_GxM_abundance_normalized_interaction_effects_PCS1and2_Plot.pdf")
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

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/IBD-eQTL_ileum_GxM_abundance_normalized_interaction_effects_taxa_labeled_PCS1and2_Plot.pdf")
pdf(fname)
pca <- ggplot(meta, aes(PC1, PC2, color=names, shape=cluster_id)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  coord_fixed()
print(pca) 
dev.off()

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/interaction_plots/IBD-eQTL_ileum_GxM_abundance_normalized_interaction_effects_taxa_labeled_no_legend_PCS1and2_Plot.pdf")
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


cv <- read_delim("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/IBD-eQTL_Ileum_covariates_with_microbial_data_10-2-2024.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cv <- as.data.frame(cv)


#subset columns for microbe abundance

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

bk1 <- c(seq(-1, -0.01, by = 0.1), -0.001)  # Negative values (down to slightly below 0)
bk2 <- c(0.001, seq(0.01, 1, by = 0.1))    # Positive values (from slightly above 0)
bk <- c(bk1, bk2)  # Combine break limits

# Create custom palette: blue for negative, gray for zero, red for positive
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1) - 1),
                "gray90", "gray90",  # Assign gray for near-zero correlations
                colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2) - 1))




png("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/GxM_microbe_abundance_correlation_heatmap.png") 
p <- pheatmap(res_sc_cor_val, color = my_palette,breaks = bk, ontsize = 9)
print(p) 
dev.off()

png("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/GxM_microbe_abundance_signficant_correlation_heatmap.png") 
p <- pheatmap(res_cor_val_sig, color = my_palette, breaks = bk, fontsize = 9)
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

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/IBD-eQTL_ileum_GxM_CLR_abundance_PCS1and2_Plot.pdf")
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

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/IBD-eQTL_ileum_GxM_CLR_abundance_taxa_labeled_PCS1and2_Plot.pdf")
pdf(fname)
pca <- ggplot(meta, aes(PC1, PC2, color=names, shape=cluster_id)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_classic() +
  coord_fixed()
print(pca) 
dev.off()

fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/plots/IBD-eQTL_ileum_GxM_CLR_abundance_taxa_labeled_no_legend_PCS1and2_Plot.pdf")
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



