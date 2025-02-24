# this script loads the eQTL mapping results for the best # of GEPCs and overlaps them with GTEx v8 results and previous eGenes
# 10/1/2019
# check for overlap between IBD rectum, ileum and GTEx colon and ileum datasets
#11 - 9 -24 SN


library(data.table)
library(qqman)
library(qvalue)
library(tidyverse)
require("VennDiagram")
library(scales)
library(UpSetR)
library(ggplot2)

FDR = 0.1


rectum_egenes <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/IBD-eQTL_rectum-permutation-egenes.txt",header=F, sep="\n", stringsAsFactors=FALSE)
rectum_egenes <- as.character(unique(rectum_egenes$V1))

rectum_all_genes <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/IBD-eQTL_rectum-all-genes.txt",header=F, sep="\n", stringsAsFactors=FALSE)
rectum_all_genes <- as.character(unique(rectum_all_genes$V1))


ileum_egenes <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/permutations/IBD-eQTL_ileum-permutation-egenes.txt",header=F, sep="\n", stringsAsFactors=FALSE)
ileum_egenes <- as.character(unique(ileum_egenes$V1))

ileum_all_genes <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/permutations/IBD-eQTL_ileum-all-genes.txt",header=F, sep="\n", stringsAsFactors=FALSE)
ileum_all_genes <- as.character(unique(ileum_all_genes$V1))


table(rectum_all_genes %in% ileum_all_genes)
# FALSE  TRUE
#   209 15179


table(rectum_egenes %in% ileum_all_genes)
# FALSE  TRUE
#    48  3407

table(ileum_egenes %in% rectum_all_genes)
# FALSE  TRUE
#    54  3258

table(ileum_egenes %in% rectum_egenes)

# FALSE  TRUE
#  1514  1798


#keep genes tested in both and egenes tested in both

allgenes <- intersect(rectum_all_genes, ileum_all_genes)
length(allgenes)
# [1] 15179


rectum_egenes <- rectum_egenes[rectum_egenes %in% allgenes]
length(rectum_egenes)
# [1] 3407

ileum_egenes <- ileum_egenes[ileum_egenes %in% allgenes]
length(ileum_egenes)
# [1] 3258


#sigmoid colon GTEx

GTEx_sigmoid <- read.table("/wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/Colon_Sigmoid.v8.egenes.txt.gz", sep="\t", header=TRUE)
GTExallgenes_sigmoid <-unique(GTEx_sigmoid$gene_id)
GTExallgenes_sigmoid <-gsub("[.].*$", "", GTExallgenes_sigmoid)
# subset to threshold:
GTEx_sigmoid_signif <- GTEx_sigmoid[GTEx_sigmoid$qval<FDR,]

GTEx_sigmoid_egenes <- as.character(unique(GTEx_sigmoid_signif$gene_id))
GTEx_sigmoid_egenes <- gsub("[.].*$", "", GTEx_sigmoid_egenes)

GTEx_sigmoid_egenesinallgenes <- GTEx_sigmoid_egenes[GTEx_sigmoid_egenes %in% allgenes]

#transverse colon GTEx
GTEx_transverse <- read.table("/wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.egenes.txt.gz", sep="\t", header=TRUE)
GTExallgenes_transverse <-unique(GTEx_transverse$gene_id)
GTExallgenes_transverse <-gsub("[.].*$", "", GTExallgenes_transverse)
# subset to threshold:
GTEx_transverse_signif <- GTEx_transverse[GTEx_transverse$qval<FDR,]

GTEx_transverse_egenes <- as.character(unique(GTEx_transverse_signif$gene_id))
GTEx_transverse_egenes <- gsub("[.].*$", "", GTEx_transverse_egenes)

GTEx_transverse_egenesinallgenes <- GTEx_transverse_egenes[GTEx_transverse_egenes %in% allgenes]

#ileum GTEx
GTEx_ileum <- read.table("/wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.egenes.txt.gz", sep="\t", header=TRUE)
GTExallgenes_ileum <-unique(GTEx_ileum$gene_id)
GTExallgenes_ileum <-gsub("[.].*$", "", GTExallgenes_ileum)
# subset to threshold:
GTEx_ileum_signif <- GTEx_ileum[GTEx_ileum$qval<FDR,]

GTEx_ileum_egenes <- as.character(unique(GTEx_ileum_signif$gene_id))
GTEx_ileum_egenes <- gsub("[.].*$", "", GTEx_ileum_egenes)

GTEx_ileum_egenesinallgenes <- GTEx_ileum_egenes[GTEx_ileum_egenes %in% allgenes]


#find egenes common to all 5 datasets
egenes_in_all <- Reduce(intersect,list(rectum_egenes, ileum_egenes, GTEx_sigmoid_egenesinallgenes, GTEx_transverse_egenesinallgenes, GTEx_ileum_egenesinallgenes))
length(egenes_in_all)
# [1] 1369


#use UpSetR to make plot to compare all 5 datasets
library(UpSetR)
library(ggplot2)

# Combine all unique genes across lists
all_genes <- unique(c(rectum_egenes, ileum_egenes, GTEx_sigmoid_egenesinallgenes, GTEx_transverse_egenesinallgenes, GTEx_ileum_egenesinallgenes))

# Create a data frame with binary indicators for each gene across all lists
gene_data <- data.frame(
  gene = all_genes,
  IBD_rectum = as.integer(all_genes %in% rectum_egenes),
  IBD_ileum = as.integer(all_genes %in% ileum_egenes),
  GTEx_sigmoid_colon = as.integer(all_genes %in% GTEx_sigmoid_egenesinallgenes),
  GTEx_transverse_colon = as.integer(all_genes %in% GTEx_transverse_egenesinallgenes),
  GTEx_terminal_ileum = as.integer(all_genes %in% GTEx_ileum_egenesinallgenes)
)

# Remove the gene column as UpSetR uses only binary columns for plotting
gene_data_upset <- gene_data[, -1]
rownames(gene_data_upset) <- gene_data$gene

colnames(gene_data_upset) <- c("IBD Rectum", "IBD Ileum", "GTEx Sigmoid Colon", "GTEx Transverse Colon", "GTEx Terminal Ileum")
mb_ratio1 <- c(0.55,0.45)
text_scale_options3 <- c(1.5, 1.25, 1.25, 1, 2, 1.5)



png("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_overlap/GTEx_IBD_eQTL_upset_plot.png", height = 2000, width = 3000, res = 300)
p <- upset(gene_data_upset, 
      sets = c("IBD Rectum", "IBD Ileum", "GTEx Sigmoid Colon", "GTEx Transverse Colon", "GTEx Terminal Ileum"),
      mb.ratio = mb_ratio1, 
	mainbar.y.label = "eGenes Shared Between Datasets", 
	 show.numbers = "yes", number.angles = 30,
	sets.x.label = "# of eGenes",
      sets.bar.color = "#00798c",
      order.by = "freq",
      main.bar.color = "#d1495b",
      matrix.color = "#0072B2",
      text.scale = text_scale_options3)
print(p)
dev.off()


#make forrest plot showing OR for each dataset enriched in their specific relevant GTEx dataset



# #   Sigmoid Colon GTEx eGene Not Sigmoid Colon GTEX eGene
# # IBD Rectum eGene                         2437                         1018
# # not IBD Rectum eGene                     5912                         6021



# fisher.test(dat)
# #  Fisher's Exact Test for Count Data

# # data:  dat
# # p-value < 2.2e-16
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  2.24579 2.64761
# # sample estimates:
# # odds ratio
# #   2.437916

sigmoid_fisher <- c("IBD Rectum In GTEx Sigmoid Colon", 2.437916, 2.24579, 2.64761)

# #   Transverse Colon GTEx eGene
# # IBD Rectum eGene                            2747
# # not IBD Rectum eGene                        6078
# #                      Not Transverse Colon GTEX eGene
# # IBD Rectum eGene                                 708
# # not IBD Rectum eGene                            5855


# fisher.test(dat)
# #  Fisher's Exact Test for Count Data

# # data:  dat
# # p-value < 2.2e-16
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  3.412989 4.095496
# # sample estimates:
# # odds ratio
# #   3.737276

transverse_fisher <- c("IBD Rectum In GTEx Transverse Colon", 3.737276, 3.412989, 4.095496)

#   Ileum GTEx eGene Not Ileum GTEX eGene
# # [1,]             1975                 1337
# # [2,]             3604                 8592



# fisher.test(dat)
# #    Fisher's Exact Test for Count Data

# # data:  dat
# # p-value < 2.2e-16
# # alternative hypothesis: true odds ratio is not equal to 1
# # 95 percent confidence interval:
# #  3.249630 3.816724
# # sample estimates:
# # odds ratio
# #   3.521289

ileum_fisher <- c("IBD Ileum In GTEx Terminal Ileum", 3.521289, 3.249630, 3.816724)


df <- rbind(ileum_fisher, sigmoid_fisher, transverse_fisher)
colnames(df) <- c("Comparison", "estimate", "ci_l", "ci_u")

df
# sigmoid_fisher     Comparison                            estimate   ci_l
# ileum_fisher      "IBD Ileum In GTEx Terminal Ileum"    "3.521289" "3.24963"
# sigmoid_fisher    "IBD Rectum In GTEx Sigmoid Colon"    "2.437916" "2.24579"
# transverse_fisher "IBD Rectum In GTEx Transverse Colon" "3.737276" "3.412989"
#                   ci_u
# ileum_fisher      "3.816724"
# sigmoid_fisher    "2.64761"
# transverse_fisher "4.095496"

df <- as.data.frame(df)
df$estimate <- as.numeric(df$estimate)
df$ci_u <- as.numeric(df$ci_u)
df$ci_l <- as.numeric(df$ci_l)

fp <- ggplot(data=df, aes(x=Comparison, y=estimate, ymin=ci_l, ymax=ci_u)) +
        geom_pointrange() +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        ylab("Enrichment (95% CI)") + theme_bw() + theme(axis.text=element_text(size=14)) + geom_errorbar(aes(ymin=ci_l, ymax=ci_u), width=0.5) +  scale_y_continuous(labels = label_number(accuracy = 0.1))

     pdf(paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_overlap/GTEx_IBD_fisher_OR_forrestplots.pdf"))
        print(fp)
    dev.off()

#enrichment of unique rectum egenes

non_rectum_genes <- unique(c(ileum_egenes, GTEx_sigmoid_egenesinallgenes, GTEx_transverse_egenesinallgenes, GTEx_ileum_egenesinallgenes))

rectum_egenes_unique <- rectum_egenes[!rectum_egenes %in% non_rectum_genes]

#rectum egene results
res <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/results/FastQTL_results_best_19GEPCs.txt",header=T, sep="\t", stringsAsFactors=FALSE)





library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(annotables)
keytypes(org.Hs.eg.db)

    original_gene_list <- res$slope

# name the vector
names(original_gene_list) <- res$pid

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df <- res %>% filter(pid %in% rectum_egenes_unique)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$slope

# Name the vector
names(genes) <- sig_genes_df$pid

# omit NA values
genes <- na.omit(genes)

go_enrich <- enrichGO(gene = names(genes),
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

#no significant results

#ileum unique gene enrichment
non_ileum_egenes <- unique(c(rectum_egenes, GTEx_sigmoid_egenesinallgenes, GTEx_transverse_egenesinallgenes, GTEx_ileum_egenesinallgenes))

ileum_egenes_unique <- ileum_egenes[!ileum_egenes %in% non_ileum_egenes]

#rectum egene results
res <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/permutations/results/FastQTL_results_best_20GEPCs.txt",header=T, sep="\t", stringsAsFactors=FALSE)





library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(annotables)
keytypes(org.Hs.eg.db)

    original_gene_list <- res$slope

# name the vector
names(original_gene_list) <- res$pid

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df <- res %>% filter(pid %in% ileum_egenes_unique)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$slope

# Name the vector
names(genes) <- sig_genes_df$pid

# omit NA values
genes <- na.omit(genes)

go_enrich <- enrichGO(gene = names(genes),
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

#no enrichment