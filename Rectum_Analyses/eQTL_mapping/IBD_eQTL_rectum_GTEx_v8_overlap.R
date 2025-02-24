# this script loads the eQTL mapping results for the best # of GEPCs and overlaps them with GTEx v8 results and previous eGenes
# 10/1/2019
# adpated for IBD eQTL project  rectum qnorm residuals 4-10-24
library(data.table)
library(qqman)
library(qvalue)
library(tidyverse)
require("VennDiagram")

FDR = 0.1

# get the number of GEPCs:
GEPCs <- read.table("eGenes-per-GEPCs.txt", sep='\t', header=T)
bestPCs <- GEPCs[GEPCs$eGenes==max(GEPCs$eGenes),"PCs"]

res <- fread(paste0("../results/FastQTL_results_best_", bestPCs, "GEPCs.txt"), sep="\t")

myegenes <- unique(res$pid[res$bqval<FDR])
allgenes <-unique(res$pid)

write.table(myegenes, file="../IBD-eQTL_rectum-permutation-egenes.txt", sep="\n", quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(allgenes, file="../IBD-eQTL_rectum-all-genes.txt", sep="\n", quote=FALSE, col.names=FALSE, row.names=FALSE)

# load GTEx v8 data:
GTEx <- read.table("/wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/Colon_Sigmoid.v8.egenes.txt.gz", sep="\t", header=TRUE)
GTExallgenes <-unique(GTEx$gene_id)
GTExallgenes <-gsub("[.].*$", "", GTExallgenes)
# subset to threshold:
GTEx_signif <- GTEx[GTEx$qval<FDR,]
sum(GTEx$qval<0.05)
# 10550
sum(GTEx$pval_beta<0.05)
# 10132
sum(GTEx$pval_nominal<0.05)
# 24483 (all)
GTExegenes <- as.character(unique(GTEx_signif$gene_id))
GTExegenes <- gsub("[.].*$", "", GTExegenes)
## GTExeqtls <- gsub("_b37","",GTEx$variant_id)
## GTExeqtls <- gsub("_",":", GTExeqtls)

GTExegenesinallgenes <- GTExegenes[GTExegenes %in% allgenes]

sum(myegenes %in% GTExegenes)/length(myegenes)
##0.7053546 70.5%
#sum(GTExeqtls %in% myeqtls)/length(GTExeqtls)
## sum( myeqtls %in% GTExeqtls)/length( myeqtls)

pdf(paste0("./plots/Venn-subsetted_egenes_sigmoid_colon_v8_", FDR, "FDR.pdf")) 
draw.pairwise.venn(area1 = length(myegenes), area2 = length(GTExegenesinallgenes), cross.area = length(myegenes[myegenes %in% GTExegenesinallgenes]), category = c(paste0("IBD-eQTL Rectum egenes at ", FDR*100,"% FDR"), paste0("Sigmoid Colon GTEx egenes at ", FDR*100, "% FDR")), lty = rep("blank", 2), fill = c("light blue", "green"),alpha = rep(0.5, 2), cat.pos = c(0, 180))
dev.off() 

myegenesGTExgenes <-myegenes[myegenes %in% GTExallgenes]

## GTExallgenes.id <- as.character(anno[anno$ens %in% GTExallgenes,"g.id"])
## write.table(GTExallgenes)

# test enrichemnt:

dat <- matrix(c(length(myegenes[myegenes %in% GTExegenesinallgenes]),
 length(GTExegenesinallgenes)-length(myegenes[myegenes %in% GTExegenesinallgenes]),
 length(myegenes[!myegenes %in% GTExegenesinallgenes]),
 length(allgenes)-length(myegenes)-length(GTExegenesinallgenes)+length(myegenes[myegenes %in% GTExegenesinallgenes])) ,2,2 )

colnames(dat) <- c("Sigmoid Colon GTEx eGene", "Not Sigmoid Colon GTEX eGene")
rownames(dat) <- c("IBD Rectum eGene", "not IBD Rectum eGene")
dat


#   Sigmoid Colon GTEx eGene Not Sigmoid Colon GTEX eGene
# IBD Rectum eGene                         2437                         1018
# not IBD Rectum eGene                     5912                         6021



fisher.test(dat)
#  Fisher's Exact Test for Count Data

# data:  dat
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.24579 2.64761
# sample estimates:
# odds ratio
#   2.437916


length(myegenes[!myegenes %in% GTExegenesinallgenes])
#1018


pdf(paste0("./plots/Venn-egenes_sigmoid_colon_in_GTEx_v8_", FDR, "FDR.pdf"))
draw.pairwise.venn(area1 = length(myegenesGTExgenes), area2 = length(GTExegenesinallgenes), cross.area = length(myegenesGTExgenes[myegenesGTExgenes %in% GTExegenesinallgenes]), category = c(paste0("IBD-eQTL Rectum egenes at ", FDR*100, "%FDR\n (tested in GTEx)"), paste0("Sigmoid Colon GTEx egenes at ", FDR*100, "% FDR\n (tested in IBD-eQTL)")), lty = rep("blank", 2), fill = c("light blue", "green"),alpha = rep(0.5, 2), cat.pos = c(0, 180))
dev.off()

## pdf("./plots/Venn-eqtls_in_GTEx10FDR.pdf")
## draw.pairwise.venn(area1 = length(myeqtlsGTExSNPs), , area2 = length(GTExeqtls), cross.area = length(myeqtlsGTExSNPs[myeqtlsGTExSNPs %in% GTExeqtls]), category = c("ALOFT eQTLs at 10% FDR\n (tested in GTEx)", "GTEx whole blood eQTLs"), lty = rep("blank", 2), fill = c("light blue", "green"),alpha = rep(0.5, 2), cat.pos = c(0, 180))
## dev.off()


#test if egenes not detected n GTEx are enriched for DEGs for microbiome from priya et al

 library(annotables)

ibd_taxa_gene_pairs <- read.table("/wsu/home/groups/piquelab/IBD_eQTL/covariates/priya_et_al_IBD_gene-taxa_results.txt", header=T, sep="\t", stringsAsFactors=FALSE)
ibd_taxa_gene_pairs$gene[ibd_taxa_gene_pairs$gene == "1-Mar"] <- "MARCH1"

anno <- grch38 %>% select(ensgene,symbol) %>% unique() %>% as.data.frame()

anno <- anno %>% filter(symbol %in% ibd_taxa_gene_pairs$gene) %>% dplyr::rename(gene = symbol)

ibd_taxa_gene_pairs <- merge(ibd_taxa_gene_pairs, anno, "gene")
dim(ibd_taxa_gene_pairs)
#1463   11

ibd_micro_genes <- unique(ibd_taxa_gene_pairs$ensgene)


egenes_noT_GTEx <- myegenes[!myegenes %in% GTExegenesinallgenes]

df <- data.frame("IBD-Microbe_Gene" = c(length(myegenesGTExgenes[myegenesGTExgenes %in% ibd_micro_genes]), length(egenes_noT_GTEx[egenes_noT_GTEx %in% ibd_micro_genes])),
				"Not_IBD-Microbe_Gene" = c(length(myegenesGTExgenes[!myegenesGTExgenes %in% ibd_micro_genes]),length(egenes_noT_GTEx[!egenes_noT_GTEx %in% ibd_micro_genes])),
				row.names = c("Rectum_eGene_in_sigmoid_GTEx", "Rectum_eGene_not_in_Sigmoid_GTEx"))

df
#                                 IBD.Microbe_Gene Not_IBD.Microbe_Gene
# Rectum_eGene_in_sigmoid_GTEx                  231                 3120
# Rectum_eGene_not_in_Sigmoid_GTEx               71                  947



fisher.test(df)

#   Fisher's Exact Test for Count Data

# data:  df
# p-value = 0.9438
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.7458113 1.3204716
# sample estimates:
# odds ratio
#  0.9875296



# load GTEx v8 data: Transverse Colon
GTEx <- read.table("/wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/Colon_Transverse.v8.egenes.txt.gz", sep="\t", header=TRUE)
GTExallgenes <-unique(GTEx$gene_id)
GTExallgenes <-gsub("[.].*$", "", GTExallgenes)
# subset to threshold:
GTEx_signif <- GTEx[GTEx$qval<FDR,]
sum(GTEx$qval<0.05)
# 11686
sum(GTEx$pval_beta<0.05)
# 10967
sum(GTEx$pval_nominal<0.05)
# 25379 (all)
GTExegenes <- as.character(unique(GTEx_signif$gene_id))
GTExegenes <- gsub("[.].*$", "", GTExegenes)
## GTExeqtls <- gsub("_b37","",GTEx$variant_id)
## GTExeqtls <- gsub("_",":", GTExeqtls)

GTExegenesinallgenes <- GTExegenes[GTExegenes %in% allgenes]

sum(myegenes %in% GTExegenes)/length(myegenes)
#0.7950594 79.5%
#sum(GTExeqtls %in% myeqtls)/length(GTExeqtls)
## sum( myeqtls %in% GTExeqtls)/length( myeqtls)

pdf(paste0("./plots/Venn-subsetted_egenes_transverse_colon_v8_", FDR, "FDR.pdf")) 
draw.pairwise.venn(area1 = length(myegenes), area2 = length(GTExegenesinallgenes), cross.area = length(myegenes[myegenes %in% GTExegenesinallgenes]), category = c(paste0("IBD-eQTL Rectum egenes at ", FDR*100,"% FDR"), paste0("Transverse Colon GTEx egenes at ", FDR*100, "% FDR")), lty = rep("blank", 2), fill = c("light blue", "green"),alpha = rep(0.5, 2), cat.pos = c(0, 180))
dev.off() 

myegenesGTExgenes <-myegenes[myegenes %in% GTExallgenes]

## GTExallgenes.id <- as.character(anno[anno$ens %in% GTExallgenes,"g.id"])
## write.table(GTExallgenes)

# test enrichemnt:
dat <- matrix(c(length(myegenes[myegenes %in% GTExegenesinallgenes]),
 length(GTExegenesinallgenes)-length(myegenes[myegenes %in% GTExegenesinallgenes]),
 length(myegenes[!myegenes %in% GTExegenesinallgenes]),
 length(allgenes)-length(myegenes)-length(GTExegenesinallgenes)+length(myegenes[myegenes %in% GTExegenesinallgenes])) ,2,2 )

colnames(dat) <- c("Transverse Colon GTEx eGene", "Not Transverse Colon GTEX eGene")
rownames(dat) <- c("IBD Rectum eGene", "not IBD Rectum eGene")
dat

#   Transverse Colon GTEx eGene
# IBD Rectum eGene                            2747
# not IBD Rectum eGene                        6078
#                      Not Transverse Colon GTEX eGene
# IBD Rectum eGene                                 708
# not IBD Rectum eGene                            5855


fisher.test(dat)
#  Fisher's Exact Test for Count Data

# data:  dat
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  3.412989 4.095496
# sample estimates:
# odds ratio
#   3.737276



length(myegenes[!myegenes %in% GTExegenesinallgenes])
#708

pdf(paste0("./plots/Venn-egenes_transverse_colon_in_GTEx_v8_", FDR, "FDR.pdf"))
draw.pairwise.venn(area1 = length(myegenesGTExgenes), area2 = length(GTExegenesinallgenes), cross.area = length(myegenesGTExgenes[myegenesGTExgenes %in% GTExegenesinallgenes]), category = c(paste0("IBD-eQTL Rectum egenes at ", FDR*100, "%FDR\n (tested in GTEx)"), paste0("Transverse Colon GTEx egenes at ", FDR*100, "% FDR\n (tested in IBD-eQTL)")), lty = rep("blank", 2), fill = c("light blue", "green"),alpha = rep(0.5, 2), cat.pos = c(0, 180))
dev.off()


egenes_noT_GTEx <- myegenes[!myegenes %in% GTExegenesinallgenes]

df <- data.frame("IBD-Microbe_Gene" = c(length(myegenesGTExgenes[myegenesGTExgenes %in% ibd_micro_genes]), length(egenes_noT_GTEx[egenes_noT_GTEx %in% ibd_micro_genes])),
				"Not_IBD-Microbe_Gene" = c(length(myegenesGTExgenes[!myegenesGTExgenes %in% ibd_micro_genes]),length(egenes_noT_GTEx[!egenes_noT_GTEx %in% ibd_micro_genes])),
				row.names = c("Rectum_eGene_in_transverse_GTEx", "Rectum_eGene_not_in_transverse_GTEx"))

df
#                                     IBD.Microbe_Gene Not_IBD.Microbe_Gene
# Rectum_eGene_in_transverse_GTEx                  231                 3174
# # Rectum_eGene_not_in_transverse_GTEx               50                  658


fisher.test(df)
#   Fisher's Exact Test for Count Data

# data:  df
# p-value = 0.806
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  0.6939681 1.3432126
# sample estimates:
# odds ratio
#  0.9577749
