require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(tidyverse)
library(stringr)
library(reshape)
library(annotables)
require(BiocParallel)

## Gene counts: this is our data for anlaysis
countFile<- "/wsu/home/groups/piquelab/IBD_eQTL/counts/IBD-eQTL_rnaseq_HTseq_gene_counts.txt"
data <- read.table(countFile, sep=" ",as.is=T,header=1,row.names=1)

#Covariates file
cv <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_correct_subject_info_based_on_genotypes_2023_01_29.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cv <- as.data.frame(cv)

#Covariates file
cv <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_updated_12_15_22.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cv <- as.data.frame(cv)
cv <- cv %>%  filter(body_site %in% c("Ileum", "Rectum"))
cv <- cv %>%  filter(TotalReads > 10000000)
#filter out individual who has invalid age (999)
cv <- cv %>% filter(!SUBJECT_ID == "C3031")


##################################################################
## annotation Table:
grch38
anno <- grch38 %>% select(ensgene,symbol,chr) %>% filter(grch38$biotype == "protein_coding") %>% unique() %>% as.data.frame()
annoY <- anno %>% filter(chr == "Y")
head(annoY)
data  <- data  %>% filter(rownames(data) %in% annoY$ensgene)

vec <- colSums(data)


cv$Y_counts <-vec[match(cv$SampleID,names(vec))]

pdf("IBD_eQTL_RNA_reads_map_toY_chr_QC.pdf")
p<-ggplot(cv, aes(x=Y_counts, fill=sex, color=sex)) +
  geom_histogram(aes(y=..density..),position="identity", alpha=0.5) +
  geom_density(alpha=0.6)+
  ggtitle("IBD eQTL RNA Sample Sex QC Density Plot") +
  labs(x = "Reads Mapped to Y Chromosome", y = "Density") +
  theme_classic()
 print(p)
dev.off()


pdf("IBD_eQTL_RNA_reads_map_toY_chr_QC_log10.pdf")
p<-ggplot(cv, aes(x=log10(Y_counts), fill=sex, color=sex)) +
  geom_histogram(aes(y=..density..),position="identity", alpha=0.5) +
  geom_density(alpha=0.6)+
  ggtitle("IBD eQTL RNA Sample Sex QC Density Plot") +
  labs(x = "Log 10(Reads Mapped to Y Chromosome)", y = "Density") +
  theme_classic()
 print(p)
dev.off()

