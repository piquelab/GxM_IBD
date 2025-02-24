library(tidyverse)
library(qvalue)
library(ggplot2)
library(data.table)
library(stringr)
library(annotables)
library(patchwork)
library(Hmisc)
library(pheatmap)

res_dichotmous <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_presence_absence/GxE_results/Ileum_GxE_presence_absence_all_results_20PCs.txt", sep="\t", stringsAsFactors=F, header=T, comment="")



nrow(res_dichotmous)
#205344


res_dichotmous$taxa_gene_pair <- paste0(res_dichotmous$var, "_", res_dichotmous$symbol)


res_abundance <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_all_results.txt", sep="\t", stringsAsFactors=F, header=T, comment="")

nrow(res_abundance)
# 423379

res_abundance$taxa_gene_pair <- paste0(res_abundance$var, "_", res_abundance$symbol)

table(unique(res_dichotmous$var) %in% unique(res_abundance$var))

# TRUE
# 62


res_abundance_sig_fdr <- res_abundance %>% filter(res_abundance$interaction_qval < 0.1)
nrow(res_abundance_sig_fdr)
#24



res_abundance_sig_pval <- res_abundance %>% filter(res_abundance$interaction_pval < 0.05)
nrow(res_abundance_sig_pval)
#20039



res_dichotmous_sig_fdr <- res_dichotmous %>% filter(res_dichotmous$interaction_qval < 0.1)
nrow(res_dichotmous_sig_fdr)
#10

res_dichotmous_sig_pval <- res_dichotmous %>% filter(res_dichotmous$interaction_pval < 0.05)
nrow(res_dichotmous_sig_pval)
# 10056

#check if any sig abundance at FDR level are significant in dichotomous GxE at pvalue < 0.05

res_abun_dichotomous <- res_dichotmous_sig_pval %>% filter(taxa_gene_pair %in% res_abundance_sig_fdr$taxa_gene_pair)
nrow(res_abun_dichotomous)
#12


res_dichotomous_abund <- res_abundance_sig_pval %>% filter(taxa_gene_pair %in% res_dichotmous_sig_fdr$taxa_gene_pair)
nrow(res_dichotomous_abund)
#6

#make directory for overlap interaction plots + data


#save results
write.table(res_abun_dichotomous, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_overlap/GxE_abundance_fdr_sig_dichotomous_pval_sig_overlap.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(res_dichotomous_abund, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_overlap/GxE_dichotomous_fdr_sig_abundance_pval_sig_overlap.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

GxM <- c(rep("Abundance GxM", 4), rep("Presence/Absence GxM", 4))
significance <- rep(c("Significant Abundance GxM", "Presence/Absence GxM p < 0.05", "Sig Presence/Absence GxM", "Abundance GxM p < 0.05"), 2)
numGxM <- c(nrow(res_abundance_sig_fdr) - nrow(res_abun_dichotomous), nrow(res_abun_dichotomous), 0, 0, 0, 0,nrow(res_dichotmous_sig_fdr) - nrow(res_dichotomous_abund) , nrow(res_dichotomous_abund))

df <- data.frame(GxM, significance, numGxM)
df$significance <- factor(df$significance, levels = c("Significant Abundance GxM", "Presence/Absence GxM p < 0.05","Sig Presence/Absence GxM", "Abundance GxM p < 0.05"))

png("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_overlap/IBD-eQTL_GxM_barplot.png")
p <- ggplot(df, aes(fill = significance, y = numGxM, x = GxM)) +
		geom_bar(position="stack", stat="identity", width = 0.5) +
     	theme_classic() +
      	scale_fill_brewer(palette = "Paired") +
      	ylab("# of GxM") +
      	coord_flip() +
      	theme(axis.title.y = element_blank())
print(p)
dev.off()

png("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_overlap/IBD-eQTL_GxM_barplot_no_legend.png", width = 600, height = 200)
p <- ggplot(df, aes(fill = significance, y = numGxM, x = GxM)) +
		geom_bar(position="stack", stat="identity", width = 0.5) +
     	theme_classic() +
      	scale_fill_brewer(palette = "Paired") +
      	ylab("# of GxM") +
      	coord_flip() +
      	theme(legend.position = "none", axis.title.y = element_blank())
print(p)
dev.off()



