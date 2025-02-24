library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(annotables)
library(tidyverse)
library(pheatmap)
library(ggtree)
library(aplot)
library(GSVA)
library(reshape2)
library(aplot)
library(RColorBrewer)
library(scales)


sc_IBD <- read_rds("/rs/rs_grp_ibdeqtl/scIBD_gex_matrix.rds")

meta <- sc_IBD[[]]

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_abundance/GxE_results/Ileum_GxM_abundance_20PCs_sig_genes.txt"
abundance_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_presence_absence/GxE_results/Ileum_GxM_presence_absence_19PCs_sig_genes.txt"
presence_absence_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)


goi_abundance <- abundance_GxM$symbol

goi_pa <- presence_absence_GxM$symbol

#ileum INTACT GxM
fname <- "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD-eQTL_Ileum_INTACT_sig_genes_with_GxM_condensed_results.txt"
INTACT_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

ileum_gxm_intact <- INTACT_GxM$symbol

goi_ileum <- unique(c(goi_abundance, goi_pa, ileum_gxm_intact))


fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/GxE_results/GxM_abundance_19PCs_sig_genes.txt"
abundance_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/GxE_results/GxM_presence_absence_19PCs_sig_genes.txt"
presence_absence_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

goi_abundance <- abundance_GxM$symbol

goi_pa <- presence_absence_GxM$symbol

fname <- "/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/INTACT_analysis/GxM_INTACT_overlap/INTACT_sig_genes_with_GxM_condensed_results.txt"
INTACT_GxM <- read.table(fname, sep="\t", header=TRUE, stringsAsFactors=FALSE)

rectum_gxm_intact <- INTACT_GxM$symbol


goi_rectum <- unique(c(goi_abundance, goi_pa, rectum_gxm_intact))



goi_all <- unique(c(goi_ileum, goi_rectum))

length(goi_all)
#81

sc_IBD_subset = sc_IBD[goi_all]
  sc_IBD_subset$minor_cluster = factor(sc_IBD_subset$minor_cluster, 
                                    levels = sc_IBD_subset$minor_cluster %>% droplevels() %>% levels %>% rev)

dp_all <- DotPlot(sc_IBD_subset, features = goi_all, group.by = "minor_cluster", cols = c("lightgray", "red"), scale = FALSE, dot.scale = 8) + RotatedAxis()
dotplot_data <- dp_all$data
dotplot_data$features.plot <- as.character(dotplot_data$features.plot)
dotplot_data$id <- as.character(dotplot_data$id)

goi_tested <- unique(dotplot_data$features.plot)
length(goi_tested)
# [1] 72

cell_order <- levels(sc_IBD_subset$minor_cluster)
scaled_avg_exp_df <- data.frame(matrix(ncol=length(cell_order), nrow=length(goi_tested)))
colnames(scaled_avg_exp_df) <- cell_order
rownames(scaled_avg_exp_df) <- goi_tested

for(i in 1:length(goi_tested)) {

   gene <- goi_tested[i]

   for(j in 1:length(cell_order)) {
      cell <- cell_order[j]

      scaled_avg_exp_df[gene, cell] <- dotplot_data[dotplot_data$features.plot == gene & dotplot_data$id == cell,]$avg.exp.scaled
   }
}
length(goi_tested)
# [1] 62

levels(sc_IBD_subset$minor_cluster)
#   [1] "Adult glia"               "Fetal glia 3"
#   [3] "Fetal glia 2"             "Fetal glia 1"
#   [5] "Differentiating glia"     "ENCC/glia progenitor"
#   [7] "Cycling ENCC/glia"        "Neuroblast"
#   [9] "Neuronal/Branch B"        "Neuronal/Branch A"
#  [11] "LEC"                      "Fetal venous capillary"
#  [13] "Fetal arterial EC"        "Cycling EC"
#  [15] "Adult arterial EC"        "Adult venous EC (C7+)"
#  [17] "Adult venous EC (SELE+)"  "Adult arterial capillary"
#  [19] "Mesothelium"              "Mesoderm 2"
#  [21] "Mesoderm 1"               "Mature pericyte"
#  [23] "Immature pericyte"        "Cycling SMC"
#  [25] "SMC 2"                    "SMC 1"
#  [27] "Myofibroblast 2"          "Myofibroblast 1"
#  [29] "Cycling stromal"          "Reticular fibroblast"
#  [31] "Inflammatory fibroblast"  "mLTo"
#  [33] "FDC"                      "Transitional stromal"
#  [35] "Stromal 3"                "Stromal 2"
#  [37] "Fetal stromal 2"          "Stromal 1"
#  [39] "Fetal stromal 1"          "Fetal enterocyte"
#  [41] "Fetal colonocyte"         "Fetal cycling TA"
#  [43] "Fetal progenitor"         "Cycling TA"
#  [45] "TA"                       "Enterocyte"
#  [47] "Adult colonocyte"         "Pediatric colonocyte"
#  [49] "DUOX2+ epithelial"        "BEST4+ epithelial"
#  [51] "Enteroendocrine"          "M-like cell"
#  [53] "Paneth"                   "Tuft"
#  [55] "Goblet"                   "Cycling plasma"
#  [57] "Cycling GC B"             "GC B"
#  [59] "Memory B"                 "Naive B"
#  [61] "Pro-B"                    "IgA-IgG- plasma"
#  [63] "IgA+IgG+ plasma"          "IgG plasma"
#  [65] "IgA plasma"               "CD56+ SELL_low NK"
#  [67] "CD56+ SELL_high NK"       "CD16+ NK"
#  [69] "NCR+ ILC3"                "NCR- ILC3"
#  [71] "CD8+ activated T"         "CD8+ MAIT"
#  [73] "CD8+ Tc17"                "CD8+ IEL"
#  [75] "CD8+ Trm"                 "CD8+ Tem"
#  [77] "CD8+ Teff"                "CD8+ Tcm"
#  [79] "CD8+ Tn"                  "CD4+ Treg"
#  [81] "CD4+ Tfh"                 "CD4+ activated T"
#  [83] "CD4+ Th17"                "CD4+ Trm"
#  [85] "CD4+ tissue-Tcm"          "CD4+ Temra"
#  [87] "CD4+ blood-Tcm"           "CD4+ Tn"
#  [89] "Megakaryocyte"            "Mast"
#  [91] "LAMP3+ DC"                "pDC"
#  [93] "cDC2"                     "cDC1"
#  [95] "Cycling macrophage"       "AREG+ macrophage"
#  [97] "LYVE1+ macrophage"        "APOE+ macrophage"
#  [99] "Inflammatory monocyte"    "Classical monocyte"
# [101] "Non-classical monocyte"


sc_IBD_subset$major_cluster = factor(sc_IBD_subset$major_cluster, 
                                    levels = sc_IBD_subset$major_cluster %>% droplevels() %>% levels %>% rev)
 levels(sc_IBD_subset$major_cluster)
[1] "Neural"      "Endothelial" "Mesenchymal" "Epithelial"  "B_Plasma"
[6] "ILC"         "CD8T"        "CD4T"        "Myeloid"


cell_type <- c(rep("Neural", 10), rep("Endothelial", 8), rep("Mesenchymal", 21), rep("Epithelial", 16), rep("B_Plasma", 10), rep("ILC", 5), rep("CD8T", 9), rep("CD4T", 9), rep("Myeloid", 13))

meta <- data.frame(cell_order, cell_type)


rownames(meta) <- meta$cell_order

meta <- meta %>% select(-c(cell_order))

meta_gene <- as.data.frame(goi_tested)

meta_gene <- meta_gene %>% mutate(body_site = ifelse(goi_tested %in% goi_ileum, ifelse(goi_tested %in% goi_rectum, "Both", "Ileum"), "Rectum"))

rownames(meta_gene) <- meta_gene$goi_tested

meta_gene <- meta_gene %>% select(-c(goi_tested))


pdf("/rs/rs_grp_ibdeqtl/figures/GxM_gene_sc_IBD_cell_atlas_only_genes_clustered_avg_exp_heatmap.pdf", height = 12, width = 14) 
p <- pheatmap(as.matrix(scaled_avg_exp_df), fontsize = 9, main = "Heatmap of Scaled Average Expression", cluster_cols=FALSE, annotation_col = meta, annotation_row = meta_gene)

print(p) 
dev.off()

# Extract the order of clusters
row_order <- p$tree_row$order



goi_ordered_scaled_avg_exp <- goi_tested[row_order]
gxm_intact <- unique(c(rectum_gxm_intact, ileum_gxm_intact))



  major_minor_mapping = sc_IBD_subset@meta.data %>%
    dplyr::select(minor_cluster, major_cluster) %>% 
    unique %>%
    mutate(group = major_cluster, p = "group")

dp <- DotPlot(sc_IBD_subset, group.by = "minor_cluster", features = goi_ordered_scaled_avg_exp,
              cols = c("lightgray", "red"), scale = FALSE, dot.scale = 8) + RotatedAxis()

# Extract feature names from the DotPlot data
dotplot_data <- dp$data
dotplot_data$features.plot <- as.character(dotplot_data$features.plot)

# Define color mapping: Red for INTACT genes, Black otherwise
gene_color_map <- setNames(
  ifelse(unique(dotplot_data$features.plot) %in% gxm_intact, "red", "black"),
  unique(dotplot_data$features.plot)
)

dp <- dp + theme(axis.text.x = element_text(color = gene_color_map, angle = 45, hjust = 1))

p2 = ggplot(major_minor_mapping, aes(x = p, y = minor_cluster, fill = major_cluster)) + 
    geom_tile(aes(fill = major_cluster)) +
    scale_fill_brewer(palette="Set3") +
    scale_x_discrete(position = "right")+ 
    xlab(NULL)+
    ylab(NULL) + 
    theme_minimal() + 
    labs(fill = "Cell type") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  dp %>%
    insert_left(p2, width = 0.03)

ggsave("/rs/rs_grp_ibdeqtl/figures/IBD_eQTL_all_GxM_scIBD__Dotplot_split_by_cell_type_ordered_by_avg_exp.pdf", width = 20, height = 22)



> goi_ileum
#  [1] "ERCC8"    "CD33"     "DDX25"    "CNN3"     "CCZ1"     "MYO5C"
#  [7] "PLAAT2"   "ARF3"     "LINS1"    "POMZP3"   "TRIM36"   "PLCD3"
# [13] "DUSP7"    "KIF27"    "GAA"      "PSMG1"    "THAP7"    "VMAC"
# [19] "SLC22A4"  "TBKBP1"   "ARHGEF38" "PNLIPRP2" "H3C6"     "MRPS17"
# [25] "FAM118A"  "MZT2A"    "DUSP14"   "CSRNP1"   "TGOLN2"
> goi_rectum
#  [1] "ZFYVE16"   "RRN3"      "C11orf21"  "TPSG1"     "HOXA7"     "SMUG1"
#  [7] "TMX4"      "GNG11"     "TUBB2B"    "CCDC82"    "SLC66A3"   "CA3"
# [13] "RPL36AL"   "ILK"       "TMEM154"   "RPS9"      "TNFRSF10C" "ABO"
# [19] "ISX"       "SNX18"     "C11orf54"  "AP2A2"     "CEACAM16"  "LPCAT4"
# [25] "TOR1AIP1"  "THOC3"     "IGFBP2"    "SEMG1"     "WDR54"     "JMJD8"
#[31] "SARDH"     "ACOT4"     "FAM167A"





#look at union of rectum and ileum degs to see if gxm genes in general are enriched for a particular cell compartment
gxm_genes <- goi_tested

FDR <- 0.1
res_ileum <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/permutations/results/FastQTL_results_best_20GEPCs.txt",header=T, sep="\t", stringsAsFactors=FALSE)
res_rectum <-  read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/results/FastQTL_results_best_19GEPCs.txt",header=T, sep="\t", stringsAsFactors=FALSE)

tested_genes <- union(res_ileum$pid, res_rectum$pid)

res_ileum_sig <- res_ileum[res_ileum$bqval<FDR,]
res_ileum_sig <- res_ileum_sig[order(res_ileum_sig$bqval),]

res_rectum_sig <- res_rectum[res_rectum$bqval<FDR,]
res_rectum_sig <- res_rectum_sig[order(res_rectum_sig$bqval),]

anno <- grch38 %>% select(ensgene,symbol) %>% filter(grch38$biotype == "protein_coding") %>% unique() %>% as.data.frame()

anno_ileum <- anno %>% filter(ensgene %in% res_ileum_sig$pid) %>% dplyr::rename(pid = ensgene)
anno_rectum <- anno %>% filter(ensgene %in% res_rectum_sig$pid) %>% dplyr::rename(pid = ensgene)

res_ileum_sig <- merge(res_ileum_sig, anno_ileum, by="pid")
res_rectum_sig <- merge(res_rectum_sig, anno_rectum, by="pid")

egenes <- union(res_ileum_sig$symbol, res_rectum_sig$symbol)

egenes_not_gxM <- egenes[!(egenes %in% gxm_genes)]

egenes <- unique(egenes)

sc_IBD_subset_egenes = sc_IBD[egenes]

sc_IBD_subset_egenes$major_cluster = factor(sc_IBD_subset_egenes$major_cluster, 
                                    levels = sc_IBD_subset_egenes$major_cluster %>% droplevels() %>% levels %>% rev)
 levels(sc_IBD_subset_egenes$major_cluster)
avg_exp <- AverageExpression(sc_IBD_subset_egenes, group.by = "major_cluster", return.seurat = TRUE)

# Get the scaled values

scaled_avg_matrix <- GetAssayData(avg_exp, slot = "data")

# Convert to a data frame
scaled_avg_df <- as.data.frame(scaled_avg_matrix)

# Add cluster information as a column
scaled_avg_df$gene <- rownames(scaled_avg_df)
head(scaled_avg_df)
scaled_avg_df <- na.omit(scaled_avg_df)

# Get the expression matrix
expr_matrix <- GetAssayData(sc_IBD_subset_egenes, slot = "data")  # Log-normalized data (default)

Idents(sc_IBD_subset_egenes) <- "major_cluster"
# Get the cluster identities
clusters <- Idents(sc_IBD_subset_egenes)

# Calculate the percentage expression
percent_exp <- sapply(unique(clusters), function(cluster) {
  cluster_cells <- WhichCells(sc_IBD_subset_egenes, idents = cluster)  # Cells in this cluster
  rowMeans(expr_matrix[, cluster_cells] > 0) * 100  # Percentage of cells with expression > 0
})

# Convert to a data frame
percent_exp_df <- as.data.frame(percent_exp)
colnames(percent_exp_df) <- unique(clusters)  # Name the columns by cluster
percent_exp_df$gene <- rownames(expr_matrix)  # Add gene names
head(percent_exp_df)
percent_exp_df <- na.omit(percent_exp_df)

percent_exp_df <- percent_exp_df[rownames(percent_exp_df) %in% rownames(scaled_avg_df),] 

expression_threshold <- 2
percent_threshold <- 10

# Check if average expression > 1 and percent expressed > 10% within each cell type
result_df <- data.frame(matrix(ncol=length(cell_order), nrow=nrow(scaled_avg_df))) 
result_df <- scaled_avg_df > expression_threshold & percent_exp_df > percent_threshold






cell_types <-  colnames(result_df)

egenes_expressed <- rownames(scaled_avg_df)

egenes_not_gxm <- egenes_expressed[!egenes_expressed %in% gxm_genes]

gxm_genes <- gxm_genes[gxm_genes %in% egenes_expressed]

express_overlap_df <- data.frame(matrix(ncol= 7, nrow=0))

colnames(express_overlap_df) <- c("Cell_type", "num_GxM_exp", "num_GxM_not_exp", "num_egene_exp", "num_egene_not_exp", "fisher_OR", "fisher_pval")

for (i in 1:length(cell_types)) {
  cell <- cell_types[i]
  data <- result_df[, c(cell)]
  


   exp_genes <- names(data[data])

   not_exp_genes <- names(data[!data])
   gxm_exp <- length(gxm_genes[gxm_genes %in% exp_genes])
   gxm_not_exp <- length(gxm_genes[gxm_genes %in% not_exp_genes])

   egene_exp <- length(egenes_not_gxm[egenes_not_gxm %in% exp_genes])
   egene_not_exp <- length(egenes_not_gxm[egenes_not_gxm %in% not_exp_genes])
  

  dat <- matrix(c(gxm_exp, egene_exp, gxm_not_exp, egene_not_exp), nrow=2)
  rownames(dat) <- c("GxM", "Not GxM")
  colnames(dat) <- c("Expressed", "Not Expressed")


    res <- fisher.test(dat)

   vec <- c(cell, gxm_exp, gxm_not_exp, egene_exp, egene_not_exp, res$estimate, res$p.value)
   express_overlap_df[nrow(express_overlap_df) + 1,] <- vec
  
}
express_overlap_df
 express_overlap_df
#      Cell_type num_GxM_exp num_GxM_not_exp num_egene_exp num_egene_not_exp
# 1       Neural           3              69            43              4308
# 2  Endothelial           1              71            36              4315
# 3  Mesenchymal           1              71            39              4312
# 4   Epithelial           1              71            40              4311
# 5     B-Plasma           2              70            35              4316
# 6          ILC           2              70            34              4317
# 7         CD8T           2              70            37              4314
# 8         CD4T           2              70            38              4313
# 9      Myeloid           2              70            43              4308
# 10        gene          72               0          4351                 0
#           fisher_OR        fisher_pval
# 1  4.35275627621794 0.0381219132830686
# 2  1.68791883958249  0.456523035782329
# 3  1.55702134837998  0.482856915558215
# 4  1.51774914755045  0.491352079163008
# 5  3.52123836201353  0.120990560075289
# 6  3.62607473633951  0.115620090622172
# 7   3.3295887591984  0.131911421283691
# 8  3.24126465709927  0.137455226506187
# 9  2.86127207559697  0.165884098977198
# 10                0                  1



sc_IBD_subset_egenes = sc_IBD[egenes]


 sc_IBD_subset_egenes$minor_cluster = factor(sc_IBD_subset_egenes$minor_cluster, 
                                    levels = sc_IBD_subset_egenes$minor_cluster %>% droplevels() %>% levels %>% rev)

sc_IBD_subset_egenes$major_cluster = factor(sc_IBD_subset_egenes$major_cluster, 
                                    levels = sc_IBD_subset_egenes$major_cluster %>% droplevels() %>% levels %>% rev)
 major_clusters <- levels(sc_IBD_subset_egenes$major_cluster)
cell_order <- levels(sc_IBD_subset_egenes$minor_cluster)
cell_type <- c(rep("Neural", 10), rep("Endothelial", 8), rep("Mesenchymal", 21), rep("Epithelial", 16), rep("B_Plasma", 10), rep("ILC", 5), rep("CD8T", 9), rep("CD4T", 9), rep("Myeloid", 13))

meta <- data.frame(cell_order, cell_type)

avg_exp <- AverageExpression(sc_IBD_subset_egenes, group.by = "minor_cluster", return.seurat = TRUE)

# Get the scaled values

scaled_avg_matrix <- GetAssayData(avg_exp, slot = "data")

# Convert to a data frame
scaled_avg_df <- as.data.frame(scaled_avg_matrix)

# Add cluster information as a column
scaled_avg_df$gene <- rownames(scaled_avg_df)
head(scaled_avg_df)
scaled_avg_df <- na.omit(scaled_avg_df)

# Get the expression matrix
expr_matrix <- GetAssayData(sc_IBD_subset_egenes, slot = "data")  # Log-normalized data (default)

Idents(sc_IBD_subset_egenes) <- "minor_cluster"
# Get the cluster identities
clusters <- Idents(sc_IBD_subset_egenes)

# Calculate the percentage expression
percent_exp <- sapply(unique(clusters), function(cluster) {
  cluster_cells <- WhichCells(sc_IBD_subset_egenes, idents = cluster)  # Cells in this cluster
  rowMeans(expr_matrix[, cluster_cells] > 0) * 100  # Percentage of cells with expression > 0
})

# Convert to a data frame
percent_exp_df <- as.data.frame(percent_exp)
colnames(percent_exp_df) <- unique(clusters)  # Name the columns by cluster
percent_exp_df$gene <- rownames(expr_matrix)  # Add gene names
head(percent_exp_df)
percent_exp_df <- na.omit(percent_exp_df)

percent_exp_df <- percent_exp_df[rownames(percent_exp_df) %in% rownames(scaled_avg_df),] 

expression_threshold <- 2
percent_threshold <- 10

# Check if average expression > 1 and percent expressed > 10% within each cell type
result_df <- data.frame(matrix(ncol=length(cell_order), nrow=nrow(scaled_avg_df))) 
result_df <- scaled_avg_df > expression_threshold & percent_exp_df > percent_threshold



egenes_expressed <- rownames(scaled_avg_df)

egenes_not_gxm <- egenes_expressed[!egenes_expressed %in% gxm_genes]

gxm_genes <- gxm_genes[gxm_genes %in% egenes_expressed]

express_overlap_df <- data.frame(matrix(ncol= 7, nrow=0))

colnames(express_overlap_df) <- c("Cell_type", "num_GxM_exp", "num_GxM_not_exp", "num_egene_exp", "num_egene_not_exp", "fisher_OR", "fisher_pval")

for (i in 1:length(cell_types)) {
  group <- major_clusters[i]
  clusters <- meta %>% filter(cell_type %in% group) %>% select(cell_order) %>% unlist()
  group_data <- result_df[, colnames(result_df) %in% clusters]
  
  # Calculate row means within the current group and store them
   expressed_vector <- rowSums(group_data)

   exp_genes <- names(expressed_vector[expressed_vector > 0])

   not_exp_genes <- names(expressed_vector[expressed_vector == 0])
   gxm_exp <- length(gxm_genes[gxm_genes %in% exp_genes])
   gxm_not_exp <- length(gxm_genes[gxm_genes %in% not_exp_genes])

   egene_exp <- length(egenes_not_gxm[egenes_not_gxm %in% exp_genes])
   egene_not_exp <- length(egenes_not_gxm[egenes_not_gxm %in% not_exp_genes])
  

  dat <- matrix(c(gxm_exp, egene_exp, gxm_not_exp, egene_not_exp), nrow=2)
  rownames(dat) <- c("GxM", "Not GxM")
  colnames(dat) <- c("Expressed", "Not Expressed")


    res <- fisher.test(dat)

   vec <- c(group, gxm_exp, gxm_not_exp, egene_exp, egene_not_exp, res$estimate, res$p.value)
   express_overlap_df[nrow(express_overlap_df) + 1,] <- vec
  
}
express_overlap_df
# Cell_type num_GxM_exp num_GxM_not_exp num_egene_exp num_egene_not_exp
# 1       Neural           2              70            76              4275
# 2  Endothelial           3              69            71              4280
# 3  Mesenchymal           3              69            85              4266
# 4   Epithelial           2              70            91              4260
# 5     B_Plasma           2              70            68              4283
# 6          ILC           2              70            43              4308
# 7         CD8T           2              70            46              4305
# 8         CD4T           3              69            42              4309
# 9      Myeloid           4              68            94              4257
# 10        <NA>           0              72             0              4351
#           fisher_OR        fisher_pval
# 1  1.60691058966438  0.364008427014723
# 2   2.6199280521927  0.118583391171805
# 3  2.18157168389572  0.171801688412233
# 4  1.33741693105337  0.663823900554737
# 5  1.79919225128314  0.316172479092572
# 6  2.86127207559697  0.165884098977198
# 7  2.67283859436028  0.183405956890497
# 8  4.45748194396053 0.0360459260503414
# 9   2.6628788491254 0.0744036961714226



express_overlap_df <- data.frame(matrix(ncol= 7, nrow=0))

colnames(express_overlap_df) <- c("Cell_type", "num_GxM_exp", "num_GxM_not_exp", "num_egene_exp", "num_egene_not_exp", "fisher_OR", "fisher_pval")

for (i in 1:ncol(result_df)) {
  cell <- colnames(result_df)[i]
  data <- result_df[, i]
  


   exp_genes <- names(data[data])

   not_exp_genes <- names(data[!data])
   gxm_exp <- length(gxm_genes[gxm_genes %in% exp_genes])
   gxm_not_exp <- length(gxm_genes[gxm_genes %in% not_exp_genes])

   egene_exp <- length(egenes_not_gxm[egenes_not_gxm %in% exp_genes])
   egene_not_exp <- length(egenes_not_gxm[egenes_not_gxm %in% not_exp_genes])
  

  dat <- matrix(c(gxm_exp, egene_exp, gxm_not_exp, egene_not_exp), nrow=2)
  rownames(dat) <- c("GxM", "Not GxM")
  colnames(dat) <- c("Expressed", "Not Expressed")


    res <- fisher.test(dat)

   vec <- c(cell, gxm_exp, gxm_not_exp, egene_exp, egene_not_exp, res$estimate, res$p.value)
   express_overlap_df[nrow(express_overlap_df) + 1,] <- vec
  
}
express_overlap_df %>% filter(fisher_pval < 0.05)

#   Cell_type num_GxM_exp num_GxM_not_exp num_egene_exp num_egene_not_exp
# 1              LEC           3              69            48              4303
# 2 CD4+ activated T           3              69            40              4311
# 3        LAMP3+ DC           3              69            47              4304
#          fisher_OR        fisher_pval
# 1 3.89531189923094 0.0493957832557892
# 2 4.68274015620293 0.0320757663366564
# 3 3.97902123009227 0.0470232876901812

express_overlap_df
#                    Cell_type num_GxM_exp num_GxM_not_exp num_egene_exp
# 1                 Adult glia           2              70            54
# 2               Fetal glia 3           1              71            37
# 3               Fetal glia 2           2              70            36
# 4               Fetal glia 1           2              70            39
# 5       Differentiating glia           1              71            34
# 6       ENCC/glia progenitor           1              71            38
# 7          Cycling ENCC/glia           1              71            39
# 8                 Neuroblast           1              71            40
# 9          Neuronal/Branch B           1              71            41
# 10         Neuronal/Branch A           1              71            41
# 11                       LEC           3              69            48
# 12    Fetal venous capillary           2              70            41
# 13         Fetal arterial EC           2              70            46
# 14                Cycling EC           2              70            44
# 15         Adult arterial EC           1              71            35
# 16     Adult venous EC (C7+)           1              71            35
# 17   Adult venous EC (SELE+)           2              70            44
# 18  Adult arterial capillary           1              71            40
# 19               Mesothelium           1              71            38
# 20                Mesoderm 2           1              71            39
# 21                Mesoderm 1           1              71            34
# 22           Mature pericyte           1              71            46
# 23         Immature pericyte           1              71            43
# 24               Cycling SMC           1              71            34
# 25                     SMC 2           1              71            36
# 26                     SMC 1           1              71            34
# 27           Myofibroblast 2           1              71            45
# 28           Myofibroblast 1           1              71            35
# 29           Cycling stromal           1              71            36
# 30      Reticular fibroblast           2              70            36
# 31   Inflammatory fibroblast           1              71            46
# 32                      mLTo           1              71            34
# 33                       FDC           1              71            30
# 34      Transitional stromal           1              71            43
# 35                 Stromal 3           2              70            37
# 36                 Stromal 2           1              71            39
# 37           Fetal stromal 2           1              71            44
# 38                 Stromal 1           1              71            43
# 39           Fetal stromal 1           1              71            36
# 40          Fetal enterocyte           1              71            41
# 41          Fetal colonocyte           1              71            37
# 42          Fetal cycling TA           1              71            37
# 43          Fetal progenitor           1              71            33
# 44                Cycling TA           2              70            48
# 45                        TA           1              71            41
# 46                Enterocyte           1              71            37
# 47          Adult colonocyte           2              70            33
# 48      Pediatric colonocyte           1              71            35
# 49         DUOX2+ epithelial           2              70            56
# 50         BEST4+ epithelial           1              71            38
# 51           Enteroendocrine           1              71            34
# 52               M-like cell           2              70            43
# 53                    Paneth           1              71            23
# 54                      Tuft           1              71            31
# 55                    Goblet           1              71            44
# 56            Cycling plasma           2              70            38
# 57              Cycling GC B           2              70            40
# 58                      GC B           2              70            42
# 59                  Memory B           1              71            40
# 60                   Naive B           1              71            36
# 61                     Pro-B           1              71            49
# 62           IgA-IgG- plasma           2              70            27
# 63           IgA+IgG+ plasma           2              70            30
# 64                IgG plasma           2              70            27
# 65                IgA plasma           2              70            27
# 66         CD56+ SELL-low NK           2              70            35
# 67        CD56+ SELL-high NK           2              70            40
# 68                  CD16+ NK           2              70            37
# 69                 NCR+ ILC3           2              70            35
# 70                 NCR- ILC3           1              71            40
# 71          CD8+ activated T           1              71            34
# 72                 CD8+ MAIT           2              70            37
# 73                 CD8+ Tc17           2              70            37
# 74                  CD8+ IEL           2              70            36
# 75                  CD8+ Trm           2              70            36
# 76                  CD8+ Tem           2              70            38
# 77                 CD8+ Teff           2              70            38
# 78                  CD8+ Tcm           2              70            38
# 79                   CD8+ Tn           2              70            35
# 80                 CD4+ Treg           2              70            37
# 81                  CD4+ Tfh           2              70            37
# 82          CD4+ activated T           3              69            40
# 83                 CD4+ Th17           2              70            38
# 84                  CD4+ Trm           2              70            37
# 85           CD4+ tissue-Tcm           2              70            36
# 86                CD4+ Temra           2              70            39
# 87            CD4+ blood-Tcm           2              70            38
# 88                   CD4+ Tn           2              70            36
# 89             Megakaryocyte           2              70            42
# 90                      Mast           2              70            41
# 91                 LAMP3+ DC           3              69            47
# 92                       pDC           2              70            42
# 93                      cDC2           1              71            34
# 94                      cDC1           1              71            32
# 95        Cycling macrophage           1              71            38
# 96          AREG+ macrophage           1              71            33
# 97         LYVE1+ macrophage           1              71            40
# 98          APOE+ macrophage           1              71            43
# 99     Inflammatory monocyte           1              71            38
# 100       Classical monocyte           1              71            36
# 101   Non-classical monocyte           1              71            42
# 102                     gene          72               0          4351
#     num_egene_not_exp        fisher_OR        fisher_pval
# 1                4297 2.27290790336772  0.231236171224951
# 2                4314 1.64193000516818  0.465444682253754
# 3                4315 3.42281438846615  0.126422116339585
# 4                4312 3.15746595966159  0.143050371366487
# 5                4317 1.78797016909961  0.438236965496539
# 6                4313 1.59835935625811  0.474221873331914
# 7                4312 1.55702134837998  0.482856915558215
# 8                4311 1.51774914755045  0.491352079163008
# 9                4310 1.48040026568904  0.499709598629065
# 10               4310 1.48040026568904  0.499709598629065
# 11               4303 3.89531189923094 0.0493957832557892
# 12               4310 3.00214353942508  0.154382464699436
# 13               4305 2.67283859436028  0.183405956890497
# 14               4307 2.79563744555557  0.171691471054204
# 15               4316 1.73653287344424   0.44745459049295
# 16               4316 1.73653287344424   0.44745459049295
# 17               4307 2.79563744555557  0.171691471054204
# 18               4311 1.51774914755045  0.491352079163008
# 19               4313 1.59835935625811  0.474221873331914
# 20               4312 1.55702134837998  0.482856915558215
# 21               4317 1.78797016909961  0.438236965496539
# 22               4305 1.31802998825887  0.539508307648559
# 23               4308 1.41091581150286  0.516020467658825
# 24               4317 1.78797016909961  0.438236965496539
# 25               4315 1.68791883958249  0.456523035782329
# 26               4317 1.78797016909961  0.438236965496539
# 27               4306 1.34762387723217  0.531806704431531
# 28               4316 1.73653287344424   0.44745459049295
# 29               4315 1.68791883958249  0.456523035782329
# 30               4315 3.42281438846615  0.126422116339585
# 31               4305 1.31802998825887  0.539508307648559
# 32               4317 1.78797016909961  0.438236965496539
# 33               4321 2.02816706938707  0.399825706100802
# 34               4308 1.41091581150286  0.516020467658825
# 35               4314  3.3295887591984  0.131911421283691
# 36               4312 1.55702134837998  0.482856915558215
# 37               4307  1.3785427845246  0.523978112401589
# 38               4308 1.41091581150286  0.516020467658825
# 39               4315 1.68791883958249  0.456523035782329
# 40               4310 1.48040026568904  0.499709598629065
# 41               4314 1.64193000516818  0.465444682253754
# 42               4314 1.64193000516818  0.465444682253754
# 43               4318 1.84266021808352  0.428867741849504
# 44               4303  2.5604358237161  0.195236474953091
# 45               4310 1.48040026568904  0.499709598629065
# 46               4314 1.64193000516818  0.465444682253754
# 47               4318  3.7365367439887  0.110314137850321
# 48               4316 1.73653287344424   0.44745459049295
# 49               4295 2.19081883244554  0.243348055528618
# 50               4313 1.59835935625811  0.474221873331914
# 51               4317 1.78797016909961  0.438236965496539
# 52               4308 2.86127207559697  0.165884098977198
# 53               4328 2.64930627889128  0.326276974847612
# 54               4320 1.96232803762237  0.409664628951609
# 55               4307  1.3785427845246  0.523978112401589
# 56               4313 3.24126465709927  0.137455226506187
# 57               4311 3.07785432434895  0.148693781054655
# 58               4309 2.93003108483879  0.160113513513782
# 59               4311 1.51774914755045  0.491352079163008
# 60               4315 1.68791883958249  0.456523035782329
# 61               4302 1.23649430845431  0.561871338300202
# 62               4324 4.57238034177209 0.0800374284613557
# 63               4321 4.11253088740751 0.0948190977258002
# 64               4324 4.57238034177209 0.0800374284613557
# 65               4324 4.57238034177209 0.0800374284613557
# 66               4316 3.52123836201353  0.120990560075289
# 67               4311 3.07785432434895  0.148693781054655
# 68               4314  3.3295887591984  0.131911421283691
# 69               4316 3.52123836201353  0.120990560075289
# 70               4311 1.51774914755045  0.491352079163008
# 71               4317 1.78797016909961  0.438236965496539
# 72               4314  3.3295887591984  0.131911421283691
# 73               4314  3.3295887591984  0.131911421283691
# 74               4315 3.42281438846615  0.126422116339585
# 75               4315 3.42281438846615  0.126422116339585
# 76               4313 3.24126465709927  0.137455226506187
# 77               4313 3.24126465709927  0.137455226506187
# 78               4313 3.24126465709927  0.137455226506187
# 79               4316 3.52123836201353  0.120990560075289
# 80               4314  3.3295887591984  0.131911421283691
# 81               4314  3.3295887591984  0.131911421283691
# 82               4311 4.68274015620293 0.0320757663366564
# 83               4313 3.24126465709927  0.137455226506187
# 84               4314  3.3295887591984  0.131911421283691
# 85               4315 3.42281438846615  0.126422116339585
# 86               4312 3.15746595966159  0.143050371366487
# 87               4313 3.24126465709927  0.137455226506187
# 88               4315 3.42281438846615  0.126422116339585
# 89               4309 2.93003108483879  0.160113513513782
# 90               4310 3.00214353942508  0.154382464699436
# 91               4304 3.97902123009227 0.0470232876901812
# 92               4309 2.93003108483879  0.160113513513782
# 93               4317 1.78797016909961  0.438236965496539
# 94               4319 1.90060768295409  0.419344461954452
# 95               4313 1.59835935625811  0.474221873331914
# 96               4318 1.84266021808352  0.428867741849504
# 97               4311 1.51774914755045  0.491352079163008
# 98               4308 1.41091581150286  0.516020467658825
# 99               4313 1.59835935625811  0.474221873331914
# 100              4315 1.68791883958249  0.456523035782329
# 101              4309 1.44483057164693  0.507931673246437
# 102                 0                0                  1
# >
