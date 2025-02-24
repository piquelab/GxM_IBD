##make heatmap of correlation of LFC zscore for all taxa with GxM and genes with GxM
library(tidyverse)
library(Hmisc)
library(pheatmap)
library(ggplot2)
library(data.table)
library(stringr)
library(annotables)
library(RColorBrewer)
library(ggraph)
library(dplyr)

platePrefix <- "microbes_with_GxE"

topDirectory <- 'out_data_'
outDir <- paste(topDirectory, platePrefix, sep='')

##
plotsDir <- paste(outDir, '/plots', sep='')

##
statsDir <- paste(outDir, '/stats', sep='')

##
dataDir <- paste(outDir, '/data_objects', sep='')

signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/signif_interactions_list_uncorrected.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("trait", "k,", "Interactions")

signif <- signif %>% filter(Interactions > 0)
dim(signif)
#40 taxa
signif <- signif %>% distinct(trait, .keep_all = TRUE)

redundant_microbes <- c("Bacteria.Firmicutes.Clostridia.Clostridiales", "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales", "Bacteria.Proteobacteria.Deltaproteobacteria", "Bacteria.Proteobacteria.Deltaproteobacteria.Desulfovibrionales.Desulfovibrionaceae")

signif <- signif %>% filter(!trait %in% redundant_microbes)

abundance_microbes <- signif$trait

signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/signif_interactions_list_uncorrected.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("trait", "k,", "Interactions")

signif <- signif %>% filter(Interactions > 0)
dim(signif)
#10
sum(signif$Interactions)
#15

table(signif$Interactions > 1)
#2 taxa with more than 1 interaction

dichotomous_microbes <- signif$trait

microbes <- c(dichotomous_microbes, abundance_microbes)

microbes <- unique(microbes)
#lists of significant genes for both models
fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxM_abundance_19PCs_sig_genes.txt"
abundance_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxM_presence_absence_19PCs_sig_genes.txt"
presence_absence_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)


goi_abundance <- abundance_GxM$pid

goi_pa <- presence_absence_GxM$pid

gxm_genes <- unique(c(goi_abundance, goi_pa))




res_list <- list()
lfc_cols <- c()

for(i in 1:length(microbes)) {
    microbe <- microbes[i]
    fname=paste("/rs/rs_grp_ibdeqtl/deseq/out_data_microbes_with_GxE/stats/", platePrefix, "_",microbe ,"_DESeq_results", ".txt", sep="")
    res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
    res$gene <- rownames(res)
    res <- res %>% filter(gene %in% gxm_genes)
    res$microbe <- microbe
    res$zscore <- res$log2FoldChange / res$lfcSE
    res <- res %>% dplyr::rename(!!paste0("taxa_",i,"_zscore") := zscore)
    res_list[[i]] <- res
    lfc_cols <- c(lfc_cols, paste0("taxa_",i,"_zscore"))
}

taxa_id_table <- cbind(microbes, lfc_cols)
write.table(taxa_id_table, file=paste("/rs/rs_grp_ibdeqtl/deseq/out_data_microbes_with_GxE/stats/", platePrefix, "_significant_taxa_id_table", ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

deseq_taxa_id_table <- taxa_id_table


microbes
# [1] "Bacteria.Firmicutes.Negativicutes.Selenomonadales.Veillonellaceae.Veillonella"
#  [2] "Bacteria.Proteobacteria.Betaproteobacteria.Burkholderiales.Alcaligenaceae.Sutterella"
#  [3] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae..Eubacterium..hallii.group"
#  [4] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Subdoligranulum"
#  [5] "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Porphyromonadaceae.Parabacteroides"
#  [6] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae.ND3007.group"
#  [7] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Fusicatenibacter"
#  [8] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
#  [9] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Erysipelotrichaceae.UCG.003"
# [10] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae..Eubacterium..eligens.group"
# [11] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae"          
# [12] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae"
# [13] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceaencultured" 
# [14] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae"          
# [15] "Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae"
# [16] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceaencultured" 
# [17] "Bacteria.Firmicutes.Clostridia.Clostridiales.Family.XIII.Family.XIII.AD3011.group"
# [18] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Ruminiclostridium.6"
# [19] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Holdemanella"
# [20] "Bacteria.Firmicutes.Negativicutes.Selenomonadales.Acidaminococcaceae.Phascolarctobacterium"
# [21] "Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae.Haemophilus"
# [22] "Bacteria.Proteobacteria.Deltaproteobacteria.Desulfovibrionales.Desulfovibrionaceae.Bilophila"
# [23] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Butyricicoccus"


microbe_short <- c("Veillonella", "Sutterella", "Eubacterium.hallii.group", "Subdoligranulum", "Parabacteroides", "Lachnospiraceae.ND3007.group", "Fusicatenibacter",
    "Lachnospira", "Erysipelotrichaceae.UCG.003", "Eubacterium.eligens.group", "Ruminococcaceae", "Erysipelotrichaceae", "Ruminococcaceaencultured", "Lachnospiraceae", "Pasteurellaceae",
    "Lachnospiraceaencultured", "Clostridiales.Family.XIII.AD3011.group", "Ruminiclostridium.6", "Holdemanella", "Phascolarctobacterium", "Haemophilus", "Bilophila", "Butyricicoccus")


##check union of DEGs
res_list <- res_list %>% reduce(full_join, by="gene")


res_sc <- res_list[,c(lfc_cols)]

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
bk1 <- c(seq(-1,-0.01,by=0.1),-0.01)
bk2 <- c(0.01,seq(0.01,1,by=0.1))
bk <- c(bk1,bk2)  #combine the break limits for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
              "gray90", "gray90",
              c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))

fname=paste0("/rs/rs_grp_ibdeqtl/deseq/out_data_microbes_with_GxE/plots/", platePrefix,"_GxM_DESeq_zscore_correlation_heatmap_Plot.pdf")
pdf(fname)          
p <- pheatmap(res_sc_cor_val, color = my_palette, fontsize = 9)
print(p) 
dev.off()

fname=paste0("/rs/rs_grp_ibdeqtl/deseq/out_data_microbes_with_GxE/plots/", platePrefix,"_GxM_DESeq_zscore_signficicant_correlations_heatmap_Plot.pdf")
pdf(fname)  
p <- pheatmap(res_cor_val_sig, color = my_palette, fontsize = 9)
print(p) 
dev.off()

deseq_lfc_zscores <- res_sc



res_list <- list()
normalized_IE_cols <- c()

for(i in 1:length(microbes)) {
    trait <- microbes[i]
    fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_", trait,".txt")
    res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
    res <- res %>% filter(pid %in% gxm_genes)
    res$microbe <- trait
    res$ie_zscore <- res$interaction_effect / res$interaction_SE
    res <- res %>% dplyr::rename(!!paste0("taxa",i,"_GxM_IE_zscore") := ie_zscore)
    res_list[[i]] <- res
    normalized_IE_cols <- c(normalized_IE_cols, paste0("taxa",i,"_GxM_IE_zscore"))

}

taxa_id_table <- cbind(microbes, normalized_IE_cols)

res_list <- res_list %>% reduce(full_join, by="pid")

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
bk1 <- c(seq(-1,-0.01,by=0.1),-0.01)
bk2 <- c(0.01,seq(0.01,1,by=0.1))
bk <- c(bk1,bk2)  #combine the break limits for purpose of graphing

my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1)-1),
              "gray90", "gray90",
              c(colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2)-1)))



pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_incl_pa_gxm_normalized_IE_heatmap.pdf") 
p <- pheatmap(res_zscore_cor_val, color = my_palette, fontsize = 9)
print(p) 
dev.off()

pdf("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/plots/interaction_plots/GxM_all_taxa_incl_pa_gxm_normalized_IE_heatmap_significant.pdf") 
p <- pheatmap(res_zscore_cor_val_sig, color = my_palette, fontsize = 9)
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
##make scatterplot of LFC zscore versus interaction effect for gene-taxa pair

#make scatterplot for sig gxm gene-taxa pairs

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxM_abundance_19PCs_sig_genes.txt"
abundance_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxM_abundance_19PCs_sig_genes_and_taxa_summary.txt"
abundance_GxM_gene_taxa <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

abundance_GxM <- merge(abundance_GxM, abundance_GxM_gene_taxa, by="symbol")

 abundance_GxM$var
#  [1] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Butyricicoccus"
#  [2] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Butyricicoccus"
#  [3] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
#  [4] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
#  [5] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceaencultured" 
#  [6] "Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae"
#  [7] "Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae.Haemophilus"
#  [8] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceaencultured" 
#  [9] "Bacteria.Firmicutes.Clostridia.Clostridiales.Family.XIII.Family.XIII.AD3011.group"
# [10] "Bacteria.Proteobacteria.Deltaproteobacteria.Desulfovibrionales.Desulfovibrionaceae.Bilophila"
# [11] "Bacteria.Firmicutes.Negativicutes.Selenomonadales.Acidaminococcaceae.Phascolarctobacterium"
# [12] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae"          
# [13] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Butyricicoccus"
# [14] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae"          
# [15] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
# [16] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Ruminiclostridium.6"
# [17] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Butyricicoccus"
# [18] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae"          
# [19] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae"
# [20] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceaencultured" 
# [21] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Holdemanella"
# [22] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Butyricicoccus"
# [23] "Bacteria.Firmicutes.Clostridia.Clostridiales.Family.XIII.Family.XIII.AD3011.group"
# [24] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Holdemanella"
# [25] "Bacteria.Proteobacteria.Gammaproteobacteria.Pasteurellales.Pasteurellaceae"


abundance_GxM$microbe <- c("Butyricicoccus", "Butyricicoccus", "Lachnospira", "Lachnospira", "Lachnospiraceaencultured", "Pasteurellaceae", "Haemophilus", "Ruminococcaceaencultured", "Clostridiales.Family.XIII.AD3011.group", "Bilophila", "Phascolarctobacterium", "Lachnospiraceae", 
    "Butyricicoccus", "Lachnospiraceae", "Lachnospira", "Ruminiclostridium.6", "Butyricicoccus", "Ruminococcaceae", "Erysipelotrichaceae", "Lachnospiraceaencultured", "Holdemanella",
    "Butyricicoccus", "Clostridiales.Family.XIII.AD3011.group", "Holdemanella", "Pasteurellaceae")

abundance_GxM$microbe_gene <- paste0(abundance_GxM$microbe,"_", abundance_GxM$pid)
fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxM_presence_absence_19PCs_sig_genes.txt"
presence_absence_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxM_presence_absence_19PCs_sig_genes_and_taxa_summary.txt"
presence_absence_GxM_gene_taxa <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

presence_absence_GxM <- merge(presence_absence_GxM, presence_absence_GxM_gene_taxa, by="symbol")
presence_absence_GxM$var

# [1] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae.ND3007.group"
#  [2] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae..Eubacterium..eligens.group"
#  [3] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Subdoligranulum"
#  [4] "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Porphyromonadaceae.Parabacteroides"
#  [5] "Bacteria.Proteobacteria.Betaproteobacteria.Burkholderiales.Alcaligenaceae.Sutterella"
#  [6] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Fusicatenibacter"
#  [7] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
#  [8] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae..Eubacterium..hallii.group"
#  [9] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Erysipelotrichaceae.UCG.003"
# [10] "Bacteria.Firmicutes.Negativicutes.Selenomonadales.Veillonellaceae.Veillonella"
# [11] "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Porphyromonadaceae.Parabacteroides"

presence_absence_GxM$microbe <- c("Lachnospiraceae.ND3007.group", "Eubacterium.eligens.group", "Subdoligranulum", "Parabacteroides", "Sutterella", "Fusicatenibacter", "Lachnospira", "Eubacterium.hallii.group", "Erysipelotrichaceae.UCG.003", "Veillonella", "Parabacteroides")

presence_absence_GxM$microbe_gene <- paste0(presence_absence_GxM$microbe, "_", presence_absence_GxM$pid)



goi_abundance <- abundance_GxM$pid

goi_pa <- presence_absence_GxM$pid

gxm_genes <- unique(c(goi_abundance, goi_pa))


res_list <- list()

for(i in 1:length(microbes)) {
    microbe <- microbes[i]
    fname=paste("/rs/rs_grp_ibdeqtl/deseq/out_data_microbes_with_GxE/stats/", platePrefix, "_",microbe ,"_DESeq_results", ".txt", sep="")
    res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
    res$gene <- rownames(res)
    res <- res %>% filter(gene %in% gxm_genes)
    res$microbe <- microbe_short[i]
    res$microbe_gene <- paste0(res$microbe, "_", res$gene)
    res$lfc_zscore <- res$log2FoldChange / res$lfcSE
   
    res_list[[i]] <- res
}

res_all_deseq <- data.frame(matrix(ncol = ncol(res_list[[1]]), nrow = 0))
        colnames(res_all_deseq) <- colnames(res_list[[1]])

        for(n in 1:length(res_list)) {
            res_all_deseq <- rbind(res_all_deseq, res_list[[n]])
        }

res_list <- list()

for(i in 1:length(microbes)) {
    trait <- microbes[i]
    fname=paste0( fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxE_abundance_19PCs_", trait,".txt"))
    res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
    res <- res %>% filter(pid %in% gxm_genes)
    res$microbe <- microbe_short[i]
    res$ie_zscore <- res$interaction_effect / res$interaction_SE
    res$microbe_effect_zscore <- res$metagene_effect / res$metagene_SE
    res$genetic_effect_zscore <- res$dosage_effect / res$dosage_SE
    res$microbe_gene <- paste0(res$microbe, "_", res$pid)
    res_list[[i]] <- res

}

res_all_gxm <- data.frame(matrix(ncol = ncol(res_list[[1]]), nrow = 0))
        colnames(res_all_gxm) <- colnames(res_all_gxm[[1]])

        for(n in 1:length(res_list)) {
            res_all_gxm <- rbind(res_all_gxm, res_list[[n]])
        }


taxa_gene_abundance <- abundance_GxM$microbe_gene

taxa_gene_pa <- presence_absence_GxM$microbe_gene

taxa_gene_pairs <- unique(c(taxa_gene_abundance, taxa_gene_pa))

res_all_gxm_sig <- res_all_gxm %>% filter(microbe_gene %in% taxa_gene_pairs)
res_all_gxm_sig <- res_all_gxm_sig %>% mutate(ie_zscore_genetic_sign = ifelse(dosage_effect < 0, ie_zscore * -1, ie_zscore))
res_all_gxm_sig <- res_all_gxm_sig %>% mutate(abs_ie_zscore = abs(ie_zscore))
res_all_gxm_sig <- res_all_gxm_sig %>% mutate(ie_effect_sign = ifelse(interaction_effect < 0, "-", "+"))
res_all_gxm_sig <- res_all_gxm_sig %>% mutate(genetic_effect_sign = ifelse(dosage_effect < 0, "-", "+"))



res_all_deseq_sig <- res_all_deseq %>% filter(microbe_gene %in% taxa_gene_pairs)
res_all_deseq_sig <- res_all_deseq_sig %>% mutate(abs_lfc_zscore = abs(lfc_zscore))

res_all_deseq_cols <- res_all_deseq_sig %>% select(microbe, microbe_gene, lfc_zscore)

res_all_gxm_cols <- res_all_gxm_sig %>% select(microbe, microbe_gene, ie_zscore_genetic_sign)


res_all <- merge(res_all_gxm_cols, res_all_deseq_cols, by="microbe_gene")

res_all <- res_all %>% select(microbe_gene, ie_zscore_genetic_sign, lfc_zscore, microbe.x) %>% dplyr::rename(microbe_short = microbe.x)


 compare_cor <- cor.test(res_all$ie_zscore_genetic_sign, res_all$lfc_zscore, method = "pearson", alternative = "two.sided")
    cor_data <- paste0("R = ", signif(compare_cor$estimate, digits = 3), ", p-value = ", signif(compare_cor$p.value, digits = 3))

    fname=paste0("/rs/rs_grp_ibdeqtl/deseq/out_data_microbes_with_GxE/plots/", platePrefix,"_GxM_IE_zscore_mult_genetic_effect_v_DESeq_zscore_scatterplot.pdf")
pdf(fname)
p <- ggplot(res_all, aes(lfc_zscore, ie_zscore_genetic_sign)) +
  geom_point() +
  ylab("Standardized GxM Interaction Effect (multipled by genetic effect sign)") +
  xlab("Standardized LFC from DESeq") +
  geom_abline(slope=1,intercept=0, linetype = 'dashed') +
  geom_smooth(data = res_all, method=lm , color="#D55E00", se=TRUE) +
  annotate(geom = "text", label = cor_data, x = 0, y = min(res_all$ie_zscore_genetic_sign) - 1, size = 7, colour = "#D55E00") +
  theme_classic() +
  theme()
print(p) 
dev.off()


