# this script is based on ALOFT/gene_counts_to_resids_251_GTEx_thresh_coom_qnorm.R
# created 10/15/2020 JR
# Modified 3-26/2024 SN


# this script takes in the full GE data object and: creates a data object subsetted to individuals present in cv file; quantile normalizes across samples and corrects for covariates

##for IBD eQTL Ileum RNAseq samples

library(tidyverse)
library(edgeR)
library(limma)
library(annotables)

countFile<- "/wsu/home/groups/piquelab/IBD_eQTL/counts/IBD-eQTL_rnaseq_HTseq_gene_counts.txt"
data1 <- read.table(countFile, sep=" ",as.is=T,header=1,row.names=1)


# load annotation to the HTSeq file:
cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_eQTL_mapping_updated_genotypePCs_2_9_23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cov <- as.data.frame(cov)
# keep only samples selected for eQTL analysis and order by dbgap:
cov <- cov %>%  filter(TotalReads > 10000000)
#filter to include only ileum samples
cov <- cov %>%  filter(body_site == "Ileum")
#filter out individual who has invalid age (999)
cv <- cov %>% filter(!SUBJECT_ID == "C3031")

##check for individuals with more than one sample
df <- cv %>% count(SUBJECT_ID)

duplicate_ind <- df %>% filter(n > 1) %>% select(SUBJECT_ID) %>% unlist()

cv1 <- cv %>% filter(!SUBJECT_ID %in% duplicate_ind)

cv2 <- cv %>% filter(SUBJECT_ID %in% duplicate_ind)

##select samples with highest proportion of clean reads
cv2_subset <- cv2 %>% select(SUBJECT_ID, SampleID, Ratio, NonDuplicated_reads) %>% arrange(SUBJECT_ID, desc(NonDuplicated_reads))

cv2_subset_filtered <- cv2_subset[!duplicated(cv2_subset$SUBJECT_ID),]

cv2 <- cv2 %>% filter(SampleID %in% cv2_subset_filtered$SampleID)
nrow(cv2)
#4
cv <- rbind(cv1, cv2)
nrow(cv)
#82

cv <- cv[order(cv$SampleID),]

# keep only samples selected for eQTL analysis in data:
data <- data1[,cv$SampleID]

grch38
anno <- grch38 %>% select(chr, biotype, start, end, strand, ensgene) %>% filter(grch38$biotype == "protein_coding") %>% as.data.frame()
x <- c(1:22)
anno <- anno[grep(paste(x, collapse="|" ), anno$chr),]
head(anno)
data  <- data  %>% filter(rownames(data) %in% anno$ensgene)

dim(data)
#19014    86

## Normalization of data
# make edgeR object
dge <- DGEList(data)
sum(colnames(dge$counts)==cv$SampleID)

#Transform counts to counts per million
cpm <- cpm(dge)
samples <- dim(data)[2]
#82
#Remove genes that are lowly expressed
table(rowSums(dge$counts==0)==samples) #Shows how many transcripts have 0 count across all samples
# FALSE  TRUE
# 17926  1088
# as did GTEx:
keep.exprs <- rowSums(cpm>=0.1)>=(0.2*samples) & rowSums(data>=6)>=(0.2*samples) #Only keep transcripts that have cpm>0.1 in at least 20% of the samples
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]
dim(dge) #protein coding genes that pass filter
#15508 82

## # Normalize data
#quantile normalize across rows - aka across samples for each gene
mat_qnorm <- apply(dge,1,function(x){qqnorm(rank(x, ties.method = "random"), plot = F)$x})
# baseline model:
design_expanded <- model.matrix(~0+ as.factor(cv$sex) + as.numeric(cv$age) + as.factor(cv$Consent) + as.numeric(cv$genPC1) + as.numeric(cv$genPC2) +as.numeric(cv$genPC3) + as.numeric(cv$Ratio))


# Fit the linear model
model <- lm(mat_qnorm ~ design_expanded, data = data.frame(mat_qnorm, cv))

# Extract residuals
residuals <- t(residuals(model))


#save the residuals:
Res <- data.frame(residuals)
colnames(Res) <- colnames(dge)



# save the residuals
write.table(Res, "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/IBD-eQTL_ileum_residuals_q_normed.txt", sep="\t", row.names=TRUE, quote=FALSE)


wxs_cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD_eQTL_WXS_metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
wxs_cov <- as.data.frame(wxs_cov)

wxs_ileum <- wxs_cov %>% filter(submitted_subject_id %in% cv$SUBJECT_ID) %>% arrange(Run)

wxs_ileum_sampleIDs <- wxs_ileum %>% select(Run) %>% unlist()

wxs_ileum_subjectIDs <- wxs_ileum %>% select(submitted_subject_id) %>% unlist()


# make a bed file to eQTL mapping: adding DNA sampleID
colnames(Res) <- cv$SUBJECT_ID
# order by DNA subject ID
Res <- Res[,wxs_ileum$submitted_subject_id]
# 2a. load the gtf file used for HTSeq count:
library(data.table)

# 2b. add #Chr start end ID
library(dplyr)
HTSeqanno <- data.frame(anno %>% group_by(ensgene) %>% summarise(minstart=min(start),maxstop=max(end), chr=first(chr), strand=first(strand)))
rownames(HTSeqanno) <- HTSeqanno$ensgene
# add TSS and stop:
HTSeqanno[HTSeqanno$strand==1,"TSS"] <- HTSeqanno[HTSeqanno$strand==1,"minstart"]
HTSeqanno[HTSeqanno$strand==1,"STOP"] <- HTSeqanno[HTSeqanno$strand==1,"maxstop"]
HTSeqanno[HTSeqanno$strand==-1,"TSS"] <- HTSeqanno[HTSeqanno$strand==-1,"minstart"]
HTSeqanno[HTSeqanno$strand==-1,"STOP"] <- HTSeqanno[HTSeqanno$strand==-1,"maxstop"]

HTSeqanno <- HTSeqanno[rownames(Res),]
dataFQTL <- cbind(HTSeqanno$chr, HTSeqanno[,c("TSS","STOP","ensgene")], Res)
colnames(dataFQTL)[1:4] <- c("#Chr", "start", "end", "ID")
# 2c. sort by start, then by chr:
dataFQTL <- dataFQTL[order(dataFQTL$start),]
dataFQTL <- dataFQTL[order(dataFQTL[,1]),]
write.table(dataFQTL, file=paste0("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis//IBD-eQTL_ileum_covariate_corrected_qnormed.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# bgzip and index
system(paste0("bgzip /rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis//IBD-eQTL_ileum_covariate_corrected_qnormed.bed"))
system(paste0("tabix -p bed /rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis//IBD-eQTL_ileum_covariate_corrected_qnormed.bed.gz"))


