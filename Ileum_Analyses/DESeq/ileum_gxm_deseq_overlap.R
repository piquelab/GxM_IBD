library(tidyverse)
library(Hmisc)
library(pheatmap)

platePrefix <- "ileum_microbes_with_GxE"

topDirectory <- 'out_data_'
outDir <- paste(topDirectory, platePrefix, sep='')

##
plotsDir <- paste(outDir, '/plots', sep='')

##
statsDir <- paste(outDir, '/stats', sep='')

##
dataDir <- paste(outDir, '/data_objects', sep='')

signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/IBD-eQTL_Ileum_GxM_abundance_interaction_summary.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("trait", "Interactions")

signif <- signif %>% filter(Interactions > 0)
signif <- signif %>% distinct(trait, .keep_all = TRUE)


redundant_microbes <- c("Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae")
signif <- signif %>% filter(!trait %in% redundant_microbes)

abundance_microbes <- signif$trait

signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_presence_absence/IBD-eQTL_Ileum_GxM_presence_absence_interaction_summary.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("trait", "Interactions")

signif <- signif %>% filter(Interactions > 0)

dichotomous_microbes <- signif$trait

microbes <- c(dichotomous_microbes, abundance_microbes)

microbes <- unique(microbes)

#lists of significant genes for both models
fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_sig_genes.txt"
abundance_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_presence_absence/GxE_results/Ileum_GxM_presence_absence_19PCs_sig_genes.txt"
presence_absence_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)


goi_abundance <- abundance_GxM$pid

goi_pa <- presence_absence_GxM$pid

gxm_genes <- unique(c(goi_abundance, goi_pa))




res_list <- list()
lfc_cols <- c()

for(i in 1:length(microbes)) {
    microbe <- microbes[i]
    fname=paste("/rs/rs_grp_ibdeqtl/deseq/out_data_ileum_microbes_with_GxE/stats/", platePrefix, "_",microbe ,"_DESeq_results", ".txt", sep="")
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
write.table(taxa_id_table, file=paste("/rs/rs_grp_ibdeqtl/deseq/out_data_ileum_microbes_with_GxE/stats/", platePrefix, "_significant_taxa_id_table", ".txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

deseq_taxa_id_table <- taxa_id_table


microbes
#  [1] "Bacteria.Verrucomicrobia"
#  [2] "Bacteria.Proteobacteria.Deltaproteobacteria"
#  [3] "Bacteria.Firmicutes.Negativicutes.Selenomonadales.Acidaminococcaceae"
#  [4] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea"
#  [5] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
#  [6] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_NK4A136_group"
#  [7] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Ruminiclostridium_6"
#  [8] "Bacteria.Verrucomicrobia.Verrucomicrobiae.Verrucomicrobiales.Verrucomicrobiaceae.Akkermansia"
#  [9] "Bacteria.Actinobacteria.Actinobacteria.Micrococcales"
# [10] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceaencultured"
# [11] "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella_9"
# [12] "Bacteria.Firmicutes.Bacilli.Bacillales.Family_XI.Gemella"
# [13] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Anaerostipes"
# [14] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Coprococcus_1"
# [15] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnoclostridium"
# [16] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_008"
# [17] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_010"
# [18] "Bacteria.Firmicutes.Clostridia.Clostridiales.Peptostreptococcaceae.Intestinibacter"
# [19] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Clostridium_innocuum_group"
# [20] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Erysipelotrichaceae_UCG_003"
# [21] "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae.Escherichia_Shigella"

microbe_short <- c("Verrucomicrobia", "Deltaproteobacteria", "Acidaminococcaceae", "Dorea", "Lachnospira", "Lachnospiraceae_NK4A136_group",
    "Ruminiclostridium_6", "Akkermansia", "Micrococcales", "Lachnospiraceaencultured", "Prevotella_9", "Gemella", "Anaerostipes" ,"Coprococcus_1",
    "Lachnoclostridium", "Lachnospiraceae_UCG_008", "Lachnospiraceae_UCG_010", "Intestinibacter", "Clostridium_innocuum_group", "Erysipelotrichaceae_UCG_003", "Escherichia_Shigella")
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
bk1 <- c(seq(-1, -0.01, by = 0.1), -0.001)  # Negative values (down to slightly below 0)
bk2 <- c(0.001, seq(0.01, 1, by = 0.1))    # Positive values (from slightly above 0)
bk <- c(bk1, bk2)  # Combine break limits

# Create custom palette: blue for negative, gray for zero, red for positive
my_palette <- c(colorRampPalette(colors = c("darkblue", "lightblue"))(n = length(bk1) - 1),
                "gray90", "gray90",  # Assign gray for near-zero correlations
                colorRampPalette(colors = c("tomato1", "darkred"))(n = length(bk2) - 1))


fname=paste0("/rs/rs_grp_ibdeqtl/deseq/out_data_ileum_microbes_with_GxE/plots/", platePrefix,"_GxM_DESeq_zscore_correlation_heatmap_Plot.pdf")
pdf(fname)          
p <- pheatmap(res_sc_cor_val, color = my_palette, breaks = bk, fontsize = 9)
print(p) 
dev.off()

fname=paste0("/rs/rs_grp_ibdeqtl/deseq/out_data_ileum_microbes_with_GxE/plots/", platePrefix,"_GxM_DESeq_zscore_signficicant_correlations_heatmap_Plot.pdf")
pdf(fname)  
p <- pheatmap(res_cor_val_sig, color = my_palette, breaks = bk, fontsize = 9)
print(p) 
dev.off()

deseq_lfc_zscores <- res_sc


#make scatterplot for sig gxm gene-taxa pairs

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_sig_genes.txt"
abundance_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_sig_genes_and_taxa_summary.txt"
abundance_GxM_gene_taxa <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

abundance_GxM <- merge(abundance_GxM, abundance_GxM_gene_taxa, by="symbol")

 abundance_GxM$var
#  [1] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnoclostridium"
#  [2] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_010"
#  [3] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Erysipelotrichaceae_UCG_003"
#  [4] "Bacteria.Actinobacteria.Actinobacteria.Micrococcales"
#  [5] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceaencultured"
#  [6] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_008"
#  [7] "Bacteria.Actinobacteria.Actinobacteria.Micrococcales"
#  [8] "Bacteria.Firmicutes.Bacilli.Bacillales.Family_XI.Gemella"
#  [9] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Coprococcus_1"
# [10] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_010"
# [11] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_010"
# [12] "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae.Escherichia_Shigella"
# [13] "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae"
# [14] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnoclostridium"
# [15] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_UCG_010"
# [16] "Bacteria.Proteobacteria.Gammaproteobacteria.Enterobacteriales.Enterobacteriaceae.Escherichia_Shigella"
# [17] "Bacteria.Bacteroidetes.Bacteroidia.Bacteroidales.Prevotellaceae.Prevotella_9"
# [18] "Bacteria.Firmicutes.Clostridia.Clostridiales.Peptostreptococcaceae.Intestinibacter"
# [19] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Coprococcus_1"
# [20] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnoclostridium"
# [21] "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales.Erysipelotrichaceae.Clostridium_innocuum_group"
# [22] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnoclostridium"
# [23] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnoclostridium"
# [24] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Anaerostipes"

abundance_GxM$microbe <- c("Lachnoclostridium", "Lachnospiraceae_UCG_010", "Erysipelotrichaceae_UCG_003", "Micrococcales", "Lachnospiraceaencultured", "Lachnospiraceae_UCG_008", "Micrococcales" ,"Gemella", "Coprococcus_1",
    "Lachnospiraceae_UCG_010", "Lachnospiraceae_UCG_010", "Escherichia_Shigella", "Enterobacteriaceae", "Lachnoclostridium", "Lachnospiraceae_UCG_010", "Escherichia_Shigella", "Prevotella_9", "Intestinibacter", "Coprococcus_1",
    "Lachnoclostridium", "Clostridium_innocuum_group", "Lachnoclostridium", "Lachnoclostridium", "Anaerostipes")

abundance_GxM$microbe_gene <- paste0(abundance_GxM$microbe,"_", abundance_GxM$pid)
fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_presence_absence/GxE_results/Ileum_GxM_presence_absence_19PCs_sig_genes.txt"
presence_absence_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_presence_absence/GxE_results/Ileum_GxM_presence_absence_20PCs_sig_genes_and_taxa_summary.txt"
presence_absence_GxM_gene_taxa <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

presence_absence_GxM <- merge(presence_absence_GxM, presence_absence_GxM_gene_taxa, by="symbol")
presence_absence_GxM$var
# [1] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Ruminiclostridium_6"
#  [2] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_NK4A136_group"
#  [3] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospira"
#  [4] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Lachnospiraceae_NK4A136_group"
#  [5] "Bacteria.Proteobacteria.Deltaproteobacteria"
#  [6] "Bacteria.Firmicutes.Clostridia.Clostridiales.Lachnospiraceae.Dorea"
#  [7] "Bacteria.Verrucomicrobia"
#  [8] "Bacteria.Verrucomicrobia.Verrucomicrobiae.Verrucomicrobiales.Verrucomicrobiaceae.Akkermansia"
#  [9] "Bacteria.Firmicutes.Negativicutes.Selenomonadales.Acidaminococcaceae"
# [10] "Bacteria.Firmicutes.Clostridia.Clostridiales.Ruminococcaceae.Ruminiclostridium_6"

presence_absence_GxM$microbe <- c("Ruminiclostridium_6", "Lachnospiraceae_NK4A136_group", "Lachnospira", "Lachnospiraceae_NK4A136_group", "Deltaproteobacteria", "Dorea", "Verrucomicrobia", "Akkermansia", "Acidaminococcaceae", "Ruminiclostridium_6")

presence_absence_GxM$microbe_gene <- paste0(presence_absence_GxM$microbe, "_", presence_absence_GxM$pid)



goi_abundance <- abundance_GxM$pid

goi_pa <- presence_absence_GxM$pid

gxm_genes <- unique(c(goi_abundance, goi_pa))


res_list <- list()

for(i in 1:length(microbes)) {
    microbe <- microbes[i]
    fname=paste("/rs/rs_grp_ibdeqtl/deseq/out_data_ileum_microbes_with_GxE/stats/", platePrefix, "_",microbe ,"_DESeq_results", ".txt", sep="")
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
    fname=paste0( fname=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/IBD-eQTL_Ileum_GxE_abundance_20PCs_", trait, ".txt"))
    res <- read.table(fname, sep="\t", stringsAsFactors=F, header=T, comment="")
    res <- res %>% filter(pid %in% gxm_genes)
    res$microbe <- microbe_short[i]
    res$ie_zscore <- res$interaction_effect / res$interaction_SE
    res$microbe_effect_zscore <- res$metagene_effect / res$metagene_SE
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

