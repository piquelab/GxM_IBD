require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
library(annotables)
library(dplyr)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(stringr)
library(reshape)
library(readr)
library(Hmisc)
library(irlba)
require(BiocParallel)

timestamp()

cores <- 8 

platePrefix <- "microbes_with_GxE"


ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## To SampleID DESeq2 in parallel, using the
## BiocParallel library
register(MulticoreParam(cores))

## Gene counts: this is our data for anlaysis
countFile<- "/wsu/home/groups/piquelab/IBD_eQTL/counts/IBD-eQTL_rnaseq_HTseq_gene_counts.txt"
data <- read.table(countFile, sep=" ",as.is=T,header=1,row.names=1)

#Covariates file
cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD-eQTL_Rectum_microbiome_all_taxa_covariates_6-7-23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cov <- as.data.frame(cov)

data <- data[,colnames(data) %in% cov$SampleID]
n.barcodes <- dim(data)[2]


##################################################################
## annotation Table:
grch38
anno <- grch38 %>% select(ensgene,symbol,chr) %>% filter(grch38$biotype == "protein_coding") %>% unique() %>% as.data.frame()
x <- c(1:22)
anno <- anno[grep(paste(x, collapse="|" ), anno$chr),]
head(anno)
data  <- data  %>% filter(rownames(data) %in% anno$ensgene)

##################################################################

## assign variables, load data, and load experiment information
topDirectory <- 'out_data_'
outDir <- paste(topDirectory, platePrefix, sep='')
system(paste("mkdir -p",outDir))
##
plotsDir <- paste(outDir, '/plots', sep='')
plotsDir <- paste(topDirectory, platePrefix, '/plots', sep='')
system(paste("mkdir -p", plotsDir))
##
statsDir <- paste(outDir, '/stats', sep='')
system(paste("mkdir -p", statsDir))
##
dataDir <- paste(outDir, '/data_objects', sep='')
system(paste("mkdir -p", dataDir))

##################################################################
## Manual conversion to R factor objects:

cov$body_site <- factor(cov$body_site)
BodySiteLevels <- levels(cov$body_site)
cat("#", BodySiteLevels, "\n")

cov$Consent <- factor(cov$Consent)
ConsentLevels <- levels(cov$Consent)
cat("#", ConsentLevels, "\n")

cov$sex <- factor(cov$sex)
SexLevels <- levels(cov$sex)
cat("#", SexLevels, "\n")

cov$race <- factor(cov$race)
RaceLevels <- levels(cov$race)
cat("#", RaceLevels, "\n")

cov$SAMPLE_SOURCE <- factor(cov$SAMPLE_SOURCE)
SAMPLE_SOURCELevels <- levels(cov$SAMPLE_SOURCE)
cat("#", SAMPLE_SOURCELevels, "\n")

cov$AFFECTION_STATUS <- factor(cov$AFFECTION_STATUS)
AFFECTION_STATUSLevels <- levels(cov$AFFECTION_STATUS)
cat("#", AFFECTION_STATUSLevels, "\n")

cov$Disease_Status <- factor(cov$Disease_Status)
Disease_StatusLevels <- levels(cov$Disease_Status)
cat("#", Disease_StatusLevels, "\n")


cov$diagnosis <- factor(cov$diagnosis)
diagnosisLevels <- levels(cov$diagnosis)
cat("#", diagnosisLevels, "\n")

cov$is_inflamed <- factor(cov$is_inflamed)
is_inflamedLevels <- levels(cov$is_inflamed)
cat("#", is_inflamedLevels, "\n")
##################################################################
## Preparing data for DEseq:
## Combine processed data into a DESeqDataSet
## & remove genes with very low coverage/expression

res_summary <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(res_summary) <- c("Microbe", "numGenes", "numDEGs")

signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_abundance/signif_interactions_list_uncorrected.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("trait", "k,", "Interactions")

signif <- signif %>% filter(Interactions > 0)
dim(signif)
#40 taxa
signif <- signif %>% distinct(trait, .keep_all = TRUE)

redundant_microbes <- c("Bacteria.Firmicutes.Clostridia.Clostridiales", "Bacteria.Firmicutes.Erysipelotrichia.Erysipelotrichales", "Bacteria.Proteobacteria.Deltaproteobacteria", "Bacteria.Proteobacteria.Deltaproteobacteria.Desulfovibrionales.Desulfovibrionaceae")

signif <- signif %>% filter(!trait %in% redundant_microbes)

abundance_microbes <- signif$trait

signif <- read.table("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/qnorm_residuals_analysis/GxM_presence_absence/signif_interactions_list_uncorrected.txt", sep="\t", stringsAsFactors=FALSE)
colnames(signif) <- c("trait", "k,", "Interactions")

signif <- signif %>% filter(Interactions > 0)
dim(signif)
#10
sum(signif$Interactions)
#15

table(signif$Interactions > 1)
#2 taxa with more than 1 interaction

dichotomous_microbes <- signif$trait

microbes <- c(dichotomous_microbes, abundance_microbes)

microbes <- unique(microbes)

table(microbes %in% colnames(cov))

model <- "Ratio + NonDuplicated_reads + age + sex + SAMPLE_SOURCE + genPC1 + genPC2 + genPC3"
   for(i in 1: length(microbes)) {
        microbe <- microbes[i]

        cv <- cov 
        cv<- (cv[which(cv$SampleID%in% colnames(data)),]) #So cov and data match
        cv <- cv[order(cv$SampleID),]
        dataX <- data[,order(colnames(data))]
        dataX <- dataX[,colnames(dataX) %in% cv$SampleID]
        colnames(dataX)==cv$SampleID


        dds <- DESeqDataSetFromMatrix(
            countData = round(dataX),
            colData = cv,
            design = as.formula(paste0("~ ",model, " + ", microbe)))
        keep <- rowSums(counts(dds)) > 0
        dds <- dds[keep,]
        colnames(dds) <- cv$SampleID
        dim(dds)

        ## Fit the model on the whole plate
        system.time(dds <- DESeq(dds,parallel=TRUE))

        ##save
        save(dds, file=paste(dataDir, '/DESeq2_', platePrefix, '_',microbe,'_','.RData', sep=''))


        res <- results(dds, parallel=TRUE)
        summary(res)
        res <- as.data.frame(res)
        fname=paste(statsDir, '/', platePrefix, "_",microbe ,"_DESeq_results", ".txt", sep="")
        write.table(res, file=fname, quote=F, sep="\t")
        res_sig <- res %>% filter(padj <0.1)
        res_sig <- res_sig[order(res_sig$padj),]
        dim(res_sig)
        fname=paste(statsDir, '/', platePrefix, "_",microbe, "_significant_DESeq_Genes", ".txt", sep="")
        write.table(res_sig, file=fname, quote=F, sep="\t")


        vec <- c(microbe, nrow(res), nrow(res_sig))
        res_summary[nrow(res_summary) + 1, ] <- vec
    }



fname=paste(statsDir, '/', platePrefix, "_summary_of_results", ".txt", sep="")
write.table(res_summary, file=fname, quote=F, sep="\t")

