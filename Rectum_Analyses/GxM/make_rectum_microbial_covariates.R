library(tidyverse)


cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_eQTL_mapping_updated_genotypePCs_2_9_23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cov <- as.data.frame(cov)
# keep only samples selected for eQTL analysis and order by dbgap:
samples_low_reads <- c("SRR8314294","SRR8314707")
cov <- cov %>%  filter(!SampleID %in% samples_low_reads)
#filter to include only ileum and rectum samples
cov <- cov %>%  filter(body_site == "Rectum")
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
cv <- rbind(cv1, cv2)


cv <- cv[order(cv$SampleID),]

wxs_cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD_eQTL_WXS_metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
wxs_cov <- as.data.frame(wxs_cov)

wxs_rectum <- wxs_cov %>% filter(submitted_subject_id %in% cv$SUBJECT_ID) %>% arrange(Run)

wxs_rectum_sampleIDs <- wxs_rectum %>% select(Run) %>% unlist()

wxs_rectum_subjectIDs <- wxs_rectum %>% select(submitted_subject_id) %>% unlist()


##order by DNA subject ID
cv <-cv[match(wxs_rectum_subjectIDs, cv$SUBJECT_ID),]

##load meta data
rectum_meta <- read.csv("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/rectum_META.csv")


#load taxa abundances
phylum_otu <- read.csv("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/rectum_p_OTU_filt_CLR.csv")
class_otu <- read.csv("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/rectum_c_OTU_filt_CLR.csv")
order_otu <- read.csv("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/rectum_o_OTU_filt_CLR.csv")
family_otu <- read.csv("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/rectum_f_OTU_filt_CLR.csv")
genus_otu <- read.csv("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/rectum_g_OTU_filt_CLR.csv")

rectum_otu_table <- rbind(phylum_otu, class_otu, order_otu, family_otu, genus_otu)
rownames(rectum_otu_table) <- rectum_otu_table[, 1]  # Set row names to the values in the first column
rectum_otu_table <- rectum_otu_table[, -1] 

colnames(rectum_otu_table) <- rectum_meta$Participant.ID

rectum_otu_table <- as.data.frame(t(rectum_otu_table))
taxa <- colnames(rectum_otu_table)

taxa <- gsub("|", ".", taxa, fixed = TRUE)
taxa <- gsub("p_", "Bacteria.", taxa, fixed = TRUE)
taxa <- gsub("c_", "", taxa, fixed = TRUE)
taxa <- gsub("o_", "", taxa, fixed = TRUE)
taxa <- gsub("f_", "", taxa, fixed = TRUE)
taxa <- gsub("g_", "", taxa, fixed = TRUE)
colnames(rectum_otu_table) <- taxa


library(tibble)
rectum_otu_table <- rownames_to_column(rectum_otu_table, var="SUBJECT_ID")
rownames(rectum_otu_table) <- NULL

rectum_otu_table <- rectum_otu_table %>% filter(SUBJECT_ID %in% wxs_rectum_subjectIDs)


dim(rectum_otu_table)
# 61 129




priya_et_al_samples <- read.table("/wsu/home/groups/piquelab/IBD_eQTL/covariates/priya_et_al_IBD_samples.txt", header=T, sep="\t", stringsAsFactors=FALSE)

priya_et_al_samples <- priya_et_al_samples %>% filter(Biopsy_location == "Rectum")
subjectIDs <- priya_et_al_samples$Participant_ID

subjectIDs[!subjectIDs %in% wxs_rectum_subjectIDs]
 # "C3031" "H4015"

subjectIDs <- subjectIDs[subjectIDs %in% wxs_rectum_subjectIDs]

subjectIDs[!subjectIDs %in% rectum_otu_table$SUBJECT_ID]
# [1] "C3005" "C3011" "C3021" "C3027" "C3035" "H4001" "M2014" "M2039" "M2047"
# [10] "P6009" "P6018" "P6024"



##get priya et al data

library(tidyverse)
library(tidyverse)
library(irlba)
require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(Hmisc)
library(RColorBrewer)
library(reshape)


cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/HMP_IBD_RNASeq_covariates_eQTL_mapping_updated_genotypePCs_2_9_23.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
cov <- as.data.frame(cov)
# keep only samples selected for eQTL analysis and order by dbgap:
samples_low_reads <- c("SRR8314294","SRR8314707")
cov <- cov %>%  filter(!SampleID %in% samples_low_reads)
#filter to include only ileum and rectum samples
cov <- cov %>%  filter(body_site == "Rectum")
#filter out individual who has invalid age (999)
cv <- cov %>% filter(!SUBJECT_ID == "C3031")

##check for individuals with more than one sample
df <- cv %>% count(SUBJECT_ID)

duplicate_ind <- df %>% filter(n > 1) %>% select(SUBJECT_ID) %>% unlist()

cv1 <- cv %>% filter(!SUBJECT_ID %in% duplicate_ind)

cv2 <- cv %>% filter(SUBJECT_ID %in% duplicate_ind)

##select samples with highest proportion of clean reads
cv2_subset <- cv2 %>% select(SUBJECT_ID, SampleID, Ratio, NonDuplicated_reads) %>% arrange(SUBJECT_ID, desc(NonDuplicated_reads))
cv2_subset
cv2_subset_filtered <- cv2_subset[!duplicated(cv2_subset$SUBJECT_ID),]
cv2_subset_filtered


cv2 <- cv2 %>% filter(SampleID %in% cv2_subset_filtered$SampleID)
dim(cv2)

cv <- rbind(cv1, cv2)
dim(cv)
cv <- cv[order(cv$SampleID),]

wxs_cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD_eQTL_WXS_metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
wxs_cov <- as.data.frame(wxs_cov)

wxs_rectum <- wxs_cov %>% filter(submitted_subject_id %in% cv$SUBJECT_ID) %>% arrange(Run)


wxs_rectum_subjectIDs <- wxs_rectum %>% select(submitted_subject_id) %>% unlist()

##order by DNA subject ID
cv <-cv[match(wxs_rectum_subjectIDs, cv$SUBJECT_ID),]

##loadsamples from priya et al
priya_et_al_samples <- read.table("/wsu/home/groups/piquelab/IBD_eQTL/covariates/priya_et_al_IBD_samples.txt", header=T, sep="\t", stringsAsFactors=FALSE)

priya_et_al_samples <- priya_et_al_samples %>% filter(Biopsy_location == "Rectum")

#get OTU table and filter by subject 
ibd_otu_table <- read.table("/wsu/home/groups/piquelab/IBD_eQTL/covariates/all_samples_clr_t_no_contam_0.001_0.1.txt", header=T, sep="\t", stringsAsFactors=FALSE)
dim(ibd_otu_table)
#78 122



subjectIDs <- priya_et_al_samples$Participant_ID



ibd_otu_table <- ibd_otu_table %>% filter(X %in% subjectIDs)

ibd_otu_table <- ibd_otu_table %>% dplyr::rename(SUBJECT_ID = X)
dim(ibd_otu_table)


 table(colnames(rectum_otu_table) %in% colnames(ibd_otu_table))

 colnames(rectum_otu_table)[!colnames(rectum_otu_table) %in% colnames(ibd_otu_table)]

# FALSE  TRUE
#    53    76


rectum_overlap <- rectum_otu_table[,colnames(rectum_otu_table) %in% colnames(ibd_otu_table)]
ibd_overlap <- ibd_otu_table[,colnames(ibd_otu_table) %in% colnames(rectum_otu_table)]


rectum_overlap <- rectum_overlap[rectum_overlap$SUBJECT_ID %in% ibd_overlap$SUBJECT_ID,]
ibd_overlap <- ibd_overlap[ibd_overlap$SUBJECT_ID %in% rectum_overlap$SUBJECT_ID,]
#merge covariate info with OTU abundances

rownames(rectum_overlap) <- rectum_overlap$SUBJECT_ID
rectum_overlap <- rectum_overlap %>% select(-c("SUBJECT_ID"))

rownames(ibd_overlap) <- ibd_overlap$SUBJECT_ID
ibd_overlap <- ibd_overlap %>% select(-c("SUBJECT_ID"))

reorder_dataframe <- function(df1, df2) {
  # Reorder columns of df1 to match df2
  df1 <- df1[, match(colnames(df2), colnames(df1))]
  # Reorder rows of df1 to match df2
  df1 <- df1[match(rownames(df2), rownames(df1)), ]
  
  return(df1)
}


ibd_overlap <- reorder_dataframe(ibd_overlap, rectum_overlap)

  cor_mat <- cor(ibd_overlap, rectum_overlap)

  cor_vector <- diag(cor_mat)

neg_cor_taxa <- names(cor_vector[cor_vector < 0])

write.table(neg_cor_taxa,  file ="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/negative_correlated_taxa.txt",sep="\t", quote=F, row.names=F,col.names=F)

png("/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/histogram_of_correlations.png")
p <- hist(cor_vector, main = "Distribution of Correlations", xlab = "Correlation Coefficient")
print(p)
dev.off()

cv_filtered <- cv %>% filter(SUBJECT_ID %in% rectum_otu_table$SUBJECT_ID)

cov_merged <- merge(cv_filtered, rectum_otu_table, by="SUBJECT_ID")
write.table(cov_merged, file="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/IBD-eQTL_rectum_covariates_with_microbial_data_4-10-2024.txt",sep="\t",  col.names=T, row.names=F, quote=F)

write.table(taxa, file ="/rs/rs_grp_ibdeqtl/FastQTL/rectum_protein_coding_residuals/rectum_microbial_covariates/1-per-var-list.txt",sep="\t", quote=F, row.names=F,col.names=F)
