library(data.table)

gts <- read.table("./PC19_eGene_genotypes.txt", sep="\t")

## add {chr, pos, variatnname} colnames from the bed file:
colnames(gts)[1:5] <- c("chr","pos","ID", "ref", "alt")

# add dbgap.IDs as colnames as in vcf file:

subject_IDs <- read.table(file="/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD-eQTL_Rectum_DNA_SubjectIDs.txt",sep="\t", stringsAsFactors=FALSE)

library(tidyverse)
subject_IDs <- subject_IDs %>% select(everything()) %>% unlist()

colnames(gts)[6:length(colnames(gts))] <- subject_IDs

# save:
write.table(gts, "IBD-eQTL_PC19_signif_eGene_genotypes.txt", sep="\t", quote=FALSE, row.names=FALSE)