# 1/25/2019 JR
##modified by SN 11-25-23
# this file uses lm to test for genotype by presence or absence of microbe for eQTLs found when correcting for 19 PCs
library(tidyverse)
library(qvalue)
library(readr)
args = commandArgs(trailingOnly=TRUE)

trait <- "Bacteria.Proteobacteria.Gammaproteobacteria.Pseudomonadales.Pseudomonadaceae.Pseudomonas"
k=0

if (length(args)>0){
    trait <- args[1]
    }

FDR <- 0.1

# load the needed files:
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

# 2. txt file with dosages:
dosbed <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/IBD-eQTL_rectum_PC19_signif_eGene_dosages.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
dosbed <- dosbed[!duplicated(dosbed[,c("ID")]),]
rownames(dosbed) <- dosbed[,3]
dos <- dosbed[,-c(1:3)]

# 3. txt file with list of testable pairs:
pairs <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/permutations/analysis/PC19_significant_topeeQTL_pairs.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)


# 4. microbial abundance values:
cv <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/dichotomous_microbes/IBD-eQTL_Rectum_covariates_dichotomous_microbe_4-24-23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cv <- as.data.frame(cv)

#subset dosages and GE to only samples we have microbial data for

dos <- dos[,colnames(dos) %in% cv$SUBJECT_ID]
GE <- GE[,colnames(GE) %in% cv$SUBJECT_ID]


cv <-cv[match(colnames(dos), cv$SUBJECT_ID),]
rownames(cv) <- cv$SUBJECT_ID

cv <- cv %>% select(trait)
colnames(cv) <- c("var")


# subset GE df to only tested genes and sort the same way as in pairs:
expression <- GE[pairs$pid,]

# duplicate dosages that are eQTLs for several genes (and order as pairs):
dosages <- dos[pairs$sid,]

# transpose expression and dosages (because lm takes columns):
expression <- t(expression)
dosages <- t(dosages)

# make a data frame that will take the results:
pairs$Intercept_pval <- NA
pairs$dosage_pval <- NA
pairs$metagene_pval <- NA
pairs$interaction_pval <- NA
pairs$Intercept_effect <- NA
pairs$dosage_effect <- NA
pairs$metagene_effect <- NA
pairs$interaction_effect <- NA
pairs$Intercept_SE <- NA
pairs$dosage_SE <- NA
pairs$metagene_SE <- NA
pairs$interaction_SE <- NA

# loop with the best number of PCs - 0:
 for (i in 1:ncol(expression)){
      model <- lm(expression[,i]~dosages[,i]*cv[,1])
      pvalues <- summary(model)$coefficients[,4]
      pairs[i,3] <- pvalues[1]
      pairs[i,4] <- pvalues[2]
      pairs[i,5] <- pvalues[3]
      pairs[i,6] <- pvalues[4]
      effects <- summary(model)$coefficients[,1]
      pairs[i,7] <- effects[1]
      pairs[i,8] <- effects[2]
      pairs[i,9] <- effects[3]
      pairs[i,10] <- effects[4]
      ses  <- summary(model)$coefficients[,2]
      pairs[i,11] <- ses[1]
      pairs[i,12] <- ses[2]
      pairs[i,13] <- ses[3]
      pairs[i,14] <- ses[4]
      }

# multiple test correct:
pairs$Intercept_qval <- qvalue(pairs$Intercept_pval)$qvalues
pairs$dosage_qval <- qvalue(pairs$dosage_pval)$qvalues
pairs$metagene_qval <- qvalue(pairs$metagene_pval)$qvalues
pairs$interaction_qval <- qvalue(pairs$interaction_pval)$qvalues

# save the results:
write.table(pairs, file=paste0("./GxE_results/GxE_presence_absence_19PCs_", trait, ".txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# report number of significant interactions:
pairs <- pairs %>% drop_na(interaction_qval)

signif <- nrow(pairs[pairs$interaction_qval<FDR,])
tab <- t(c(trait,k, signif))

write.table(tab, file="./signif_interactions_list_uncorrected.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)

## make a qqplot and p-value histogram:
pdf(paste0("./plots/QQplots/GxE_presence_absence_19PCs_pvalues_QQplot_", trait, ".pdf"))
library(qqman)
qq(pairs$interaction_pval)
dev.off()
pdf(paste0("./plots/pval_histograms/GxE_presence_absence_19PCs_pvalues_histogram_", trait, ".pdf"))
hist(pairs$interaction_pval)
dev.off()
