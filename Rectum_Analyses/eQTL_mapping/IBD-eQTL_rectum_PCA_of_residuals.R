## this script calculates PCs for the bed file used for FastQTL analysis
#IBD eQTL Project - Rectum samples - protein coding genes only, the residuals
#11-18-23

library(tidyverse)
library(irlba)
## load normalized data / or residuals:
Res <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/IBD-eQTL_rectum_residuals_qnorm.txt", sep="\t", header=TRUE)


#make colnames the DNA subject IDs
cvf <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_eQTL_mapping_updated_genotypePCs_2_9_23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cvf <- as.data.frame(cvf)

cvf <- cvf %>% filter(SampleID %in% colnames(Res))
Res <- Res[,cvf$SampleID]

identical(cvf$SampleID, colnames(Res))

colnames(Res) <- cvf$SUBJECT_ID

wxs_cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD_eQTL_WXS_metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
wxs_cov <- as.data.frame(wxs_cov)

wxs_rectum <- wxs_cov %>% filter(submitted_subject_id %in% cvf$SUBJECT_ID) %>% arrange(Run)

Res <- Res[,wxs_rectum$submitted_subject_id]


# run PCA on residuals
PCs <- prcomp_irlba(t(Res), n=30)
summary(PCs)
mypcs <- as.data.frame(PCs$x)
rownames(mypcs) <-colnames(Res)
samples <- colnames(Res)

covs <- as_tibble(t(mypcs))
cvs <- cbind(id=colnames(mypcs),covs)



write.table(cvs,  "./IBD_eQTL_Rectum-PCcovariates-FastQTL.txt", row.names=FALSE, sep="\t", quote=FALSE)
write.table(cvs,  "./IBD_eQTL_Rectum-PCcovariates-FastQTL_nohead.txt", row.names=FALSE, sep="\t", quote=FALSE, col.names=FALSE)

PCvar <- data.frame(summary(PCs)$importance)[2,]
PCvar <- t(PCvar)

# save the proportion of variance expplained by PCs:
write.table(PCvar, "IBD-eQTL_rectum_var_explained_by_GE_PCs.txt", quote=FALSE, sep="\t")

write.table(samples, "./sample-list.txt", sep="\n", quote=FALSE, row.names=FALSE, col.names=FALSE)
