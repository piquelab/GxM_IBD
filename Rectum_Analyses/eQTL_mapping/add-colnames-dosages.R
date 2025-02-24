library(data.table)

dos <- read.table("./PC19_signif_dosages.txt", sep="\t")

## add {chr, pos, variatnname} colnames from the bed file:
colnames(dos)[1:3] <- c("chr","pos","ID")

# add dbgap.IDs as colnames as in vcf file:
# load annotation to the HTSeq file:


subject_IDs <- read.table(file="/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD-eQTL_Rectum_DNA_SubjectIDs.txt",sep="\t", stringsAsFactors=FALSE)

library(tidyverse)
subject_IDs <- subject_IDs %>% select(everything()) %>% unlist()

colnames(dos)[4:length(colnames(dos))] <- subject_IDs

# save:
write.table(dos, "IBD-eQTL_rectum_PC19_signif_eGene_dosages.txt", sep="\t", quote=FALSE, row.names=FALSE)

# dim(dos)
# [1] 3424   89

