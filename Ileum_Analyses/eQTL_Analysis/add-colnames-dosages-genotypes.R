library(data.table)

gts <- read.table("./PC20_signif_genotypes.txt", sep="\t")

## add {chr, pos, variatnname} colnames from the bed file:
colnames(gts)[1:5] <- c("chr","pos","ID", "ref", "alt")

# add subject IDs as colnames as in vcf file:

subject_IDs <- read.table(file="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/IBD-eQTL_ileum_DNA_SubjectIDs.txt",sep="\t", stringsAsFactors=FALSE)

library(tidyverse)
subject_IDs <- subject_IDs %>% select(everything()) %>% unlist()

colnames(gts)[6:length(colnames(gts))] <- subject_IDs

# save:
write.table(gts, "IBD-eQTL_Ileum_PC20_signif_eGene_genotypes.txt", sep="\t", quote=FALSE, row.names=FALSE)

library(data.table)

dos <- read.table("./PC20_signif_dosages.txt", sep="\t")

## add {chr, pos, variatnname} colnames from the bed file:
colnames(dos)[1:3] <- c("chr","pos","ID")

# add subject IDs as colnames as in vcf file:
# load annotation to the HTSeq file:

subject_IDs <- read.table(file="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/IBD-eQTL_ileum_DNA_SubjectIDs.txt",sep="\t", stringsAsFactors=FALSE)

library(tidyverse)
subject_IDs <- subject_IDs %>% select(everything()) %>% unlist()

colnames(dos)[4:length(colnames(dos))] <- subject_IDs

# save:
write.table(dos, "IBD-eQTL_ileum_PC20_signif_dosages.txt", sep="\t", quote=FALSE, row.names=FALSE)