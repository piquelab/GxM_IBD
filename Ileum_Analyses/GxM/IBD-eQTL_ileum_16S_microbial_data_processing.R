 library(tidyverse)
 library(data.table) ## for fread()
library(reshape2)
library(ggplot2)
library(RColorBrewer)


##################### Functions ####################
normalize_taxaname <- function(taxa){
  ##debug
  # microbe_name <- "Bacteria; __Firmicutes; __Clostridia; __Clostridiales; __Clostridiales_vadinBB60_group; __g" ## for debugging
  ##debug
  microbe_name <- taxa
  split <- strsplit(microbe_name,split="\\; __")
  microbe_name <- paste(split[[1]],collapse=";") # Create collapsed names
  microbe_name <- gsub("(;[a-z])+$","",microbe_name,perl=T) ## get rid of any lingering ;X
  # microbe_name <- gsub("[[:digit:]]+","",microbe_name,perl=T) ## get rid of number added at end by rowsum in summarize_taxa() to create unique groups. 
  ## These are anyways unique taxa, so we can get rid of them during visualization after downstream analysis. 
  return(microbe_name)
}

## Normalize taxa names
## Convert taxa names to form kingdom;phylum;...;species
## Addiitonal notes for IBD taxa nomenclature (SILVA DB)
## 1. For IBD (SILVA db), taxa names only go till genus level (i.e. rank = 6)
##    In this dataset, species are combined with genus label via a seperator _, 
##    instead of the separator "; __" in IBD. 
##    We don't need to sperate genus and species, rather leave them as a one entity at genus level. 
##    May be remove the separator "_" between such genera names. 
## 2. Also, some genus + species names begin with 3 _'s, e.g. "; ___Eubacterium_eligens_group".
##    This indicates that the genus is uncharacterized, so should be [Eubacterium] eligens group
##    13 such taxa: grep(";_", rownames(taxa), value = T) 
normalize_taxaname_IBD <- function(taxa){
  ##debug
  # microbe_name <- "Bacteria; __Firmicutes; __Clostridia; __Clostridiales; __Family_XIII; ___Eubacterium_brachy_group" ## for debugging
  ## debug
  microbe_name <- taxa
  split <- strsplit(microbe_name,split="\\; __")
  microbe_name <- paste(split[[1]],collapse=";") # Create collapsed names
  microbe_name <- gsub("(;[a-z])?","",microbe_name,perl=T) ## get rid of any lingering ;o;f;g ,etc.
  # microbe_name <- gsub("[A-z]__","",microbe_name,perl=T) ## get rid of initial k__ -- not applicable to IBD
  
  ## Address issue 2, i.e. genus + species names begin with 3 _'s, e.g. "; ___Eubacterium_eligens_group".
  if(grepl(";_",microbe_name)){ 
    ## I am assuming this happens at genus level
    genus_name <- tail(strsplit(microbe_name,split="\\;")[[1]],1)
    genus_name <- gsub("^_","",genus_name)
    genus_name_split <- strsplit(genus_name, split="_")
    genus_name_split[[1]][1] <- paste0("[",genus_name_split[[1]][1],"]")
    genus_name <- paste(genus_name_split[[1]], collapse = " ")
    microbe_name <- paste(head(strsplit(microbe_name,split="\\;")[[1]],5),collapse = ";")
    microbe_name <- paste(microbe_name, genus_name,sep =";")
  }
  
  ## Get rid of lingering "_" in the taxaname overall 
  microbe_name <- gsub("_"," ",microbe_name)
  
  return(microbe_name)
}

## clr + imputation function
## Gabe's function (from Dan Knight's lab)
clr_transform <- function(taxa){
  clr.taxa <- taxa
  clr.taxa = t(clr.taxa); eps = 0.5 ## changing eps alters the number of columns filterd
  clr.taxa = clr.taxa*(1 - rowSums(clr.taxa==0)*eps/rowSums(clr.taxa))
  clr.taxa[clr.taxa==0]=eps
  clr.taxa = sweep(clr.taxa,1,rowSums(clr.taxa),'/');
  ls = log(clr.taxa)
  clr.taxa = t(ls - rowMeans(ls))
  clr.taxa = clr.taxa[,!is.nan(colSums(clr.taxa))]
  return(clr.taxa)
}

## collapse at different taxa levels
summarize_taxa_IBD <- function(otu,num_samples,sep){
  lsTaxa <- list()
  for(rank in 2:6){ ## IBD data doesn't have 7th (i.e. species level)
    
    
    ## Split taxa names by seperator 
    split = strsplit(as.character(otu$taxonomy),sep) ## list of lists
    ## generate taxa name at a given taxa rank
    taxaStrings = sapply(split,function(x) paste(x[1:rank],collapse=";"))
    ## get rid of k__ at begining of string
    #taxaStrings <- gsub('([[:alpha:]])__','', taxaStrings)
    ## Remove any trailing ;
    #taxaStrings = gsub(";+$","",taxaStrings,perl=T)# Clean tips
    ## Collapse table grouped by taxa names
    taxa = rowsum(otu[,-ncol(otu)], group =  taxaStrings) 
    dim(taxa)
    
    ## remove taxa with trailing NAs
    select <- grep("(;NA)+$",rownames(taxa))
    if(length(select)!= 0){
      taxa <- taxa[-select,]
    }
    dim(taxa)
    
    ## Prevalence-based filtering
    ## Tune the prevalence filter to find out which one works for this data. 
    # select <- rowSums(taxa > 0) >= num_samples*0.5 #Keep otus which are present in more than 50% of the samples.
    # taxa <- taxa[select,]
    
    ## Filter further by only keeping taxa with rel. abundance atleast 0.1% (or 0.001) in atleast 10% of the samples
    taxa.rel <- sweep(taxa,2,colSums(taxa),'/')
    select <- rowSums(taxa.rel >= 0.001) >= num_samples*0.1
    # select <- rowSums(taxa.rel >= 0.0001) >= num_samples*0.1
    taxa <- taxa[select,]
    dim(taxa)
    
    lsTaxa[[rank]] <- taxa
  }
  return(lsTaxa)
}


sample.info <- read.csv(file="/rs/rs_grp_ibdeqtl/microbial_data_processed_by_sabrina/hmp2_metadata_2018-08-20.csv")
dim(sample.info)
table(sample.info$biopsy_location)

#otu
otu <- read.csv(file="/rs/rs_grp_ibdeqtl/microbial_data_processed_by_sabrina/otu.csv")
taxa <- read.csv(file="/rs/rs_grp_ibdeqtl/microbial_data_processed_by_sabrina/taxa.csv")

otu <- merge(otu, taxa, by = "X") %>%
column_to_rownames("X") %>%
  dplyr::rename("taxonomy"=Rank1)

  otus.to.filter <- rownames(otu[grep("Archaea", otu$taxonomy),])

  otu <- otu[!(rownames(otu) %in% otus.to.filter),]
dim(otu) 


## normalize taxa labels to trim-off stuff at uncharacterized levels 
otu.labels.fixed <- as.vector(sapply(as.character(otu$taxonomy),normalize_taxaname_IBD))

## For IBD, add the normalized taxanames to the OTU table
otu$taxonomy <- otu.labels.fixed


fileName <- "/rs/rs_grp_ibdeqtl/microbial_data_processed_by_sabrina/contaminants_all_taxa_levels_V5.txt"
contaminants <- readLines(paste0(fileName))
## replace special characters like [,] with \\ to allow grep to work appropriately
contaminants <- gsub("\\[","\\\\[",contaminants, perl = T)
contaminants <- gsub("\\]","\\\\]",contaminants, perl = T)
## extract star-marked contaminant from list (flavobacteriia and sphingobacteriia need special treatment)
special_contam <- contaminants[grep("^\\*", contaminants)]
contaminants <- contaminants[-(grep("^\\*", contaminants))]
# grep("Flavobacteriia",rownames(otu), value = T)
## filter out all contam except flavoobacteriia lineages
otu.filt <- otu
dim(otu.filt) 
for(i in 1:length(contaminants)){
  # i <- 2
  select <- grep(contaminants[i], otu.filt$taxonomy) ## This assumes that these labels are in same order as orginally in OTU table 
  if(length(select) != 0){
    otu.filt <- otu.filt[-select,]
    print(paste0("Removed ",contaminants[i],":",dim(otu.filt)[1]));flush.console()
  }
}
dim(otu.filt) #[1]  820 180
## Filter out all Falvobacteriia lineages except Capnocytophaga and Elizabethkingia
select <- grepl("Flavobacteriia", otu.filt$taxonomy) & !grepl("Capnocytophaga",otu.filt$taxonomy) & !grepl("Elizabethkingia",otu.filt$taxonomy)
length(which(select)) #13
otu.filt <- otu.filt[!select,]
grep("Flavobacteriia",otu.filt$taxonomy, value = T) ##only 3 taxa matching 
dim(otu.filt) 
# 807 180

## Filter out all Sphingobacteriia lineages except Sphingobacterium
select <- grepl("Sphingobacteriia",otu.filt$taxonomy) & !grepl("Sphingobacterium",otu.filt$taxonomy)
length(which(select)) #19. includes Pedobacter, Parapedobacter, Olivibacter, etc. 
otu.filt <- otu.filt[!select,]
grep("Sphingobacteriia",otu.filt$taxonomy, value = T) ##taxa from Sphingobacterium lineage
dim(otu.filt) 
 # 788 180
## Total OTUs removed from IBD.
dim(otu)[1] - dim(otu.filt)[1]
# 188

ileum.samples <- sample.info[sample.info$biopsy_location=="Ileum" & sample.info$data_type == "biopsy_16S",]$External.ID
length(ileum.samples) 

ileum_meta <- sample.info %>% filter(sample.info$biopsy_location=="Ileum" & sample.info$data_type == "biopsy_16S")

#filter taxa information to have just ileum samples of interest

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

wxs_cov <- read_delim("/wsu/home/groups/piquelab/IBD_eQTL/covariates/IBD_eQTL_WXS_metadata.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
wxs_cov <- as.data.frame(wxs_cov)

wxs_ileum <- wxs_cov %>% filter(submitted_subject_id %in% cv$SUBJECT_ID) %>% arrange(Run)


ileum_meta <- ileum_meta %>% filter(Participant.ID %in% cv$SUBJECT_ID)

colnames(otu.filt) <- sub("X", "", colnames(otu.filt))

cols <- c(ileum_meta$External.ID, "taxonomy")

otu.filt <- otu.filt[,cols]

dim(otu.filt)
# [1] 788  73



#taxa.levels <- summarize_taxa_IBD(otu.filt,dim(otu.filt[,-ncol(otu.filt)])[1],";") ## set taxa table specific separator

#manually do this

otu <- otu.filt
num_samples <- nrow(ileum_meta)
sep <- ";"

lsTaxa <- list()
  for(rank in 2:6){ ## IBD data doesn't have 7th (i.e. species level)
    
    
    ## Split taxa names by seperator 
    split = strsplit(as.character(otu$taxonomy),sep) ## list of lists
    ## generate taxa name at a given taxa rank
    taxaStrings = sapply(split,function(x) paste(x[1:rank],collapse=";"))
    ## get rid of k__ at begining of string
    #taxaStrings <- gsub('([[:alpha:]])__','', taxaStrings)
    ## Remove any trailing ;
    #taxaStrings = gsub(";+$","",taxaStrings,perl=T)# Clean tips
    ## Collapse table grouped by taxa names
    taxa = rowsum(otu[,-ncol(otu)], group =  taxaStrings) 
    dim(taxa)
    
    ## remove taxa with trailing NAs
    select <- grep("(;NA)+$",rownames(taxa))
    if(length(select)!= 0){
      taxa <- taxa[-select,]
    }
    dim(taxa)
    
    ## Prevalence-based filtering
    ## Tune the prevalence filter to find out which one works for this data. 
    # select <- rowSums(taxa > 0) >= num_samples*0.5 #Keep otus which are present in more than 50% of the samples.
    # taxa <- taxa[select,]
    
    ## Filter further by only keeping taxa with rel. abundance atleast 0.1% (or 0.001) in atleast 10% of the samples
    taxa.rel <- sweep(taxa,2,colSums(taxa),'/')
    select <- rowSums(taxa.rel >= 0.001) >= num_samples*0.1
    # select <- rowSums(taxa.rel >= 0.0001) >= num_samples*0.1
    taxa <- taxa[select,]
    dim(taxa)
    
    lsTaxa[[rank]] <- taxa
  }

## combine all taxa levels
#taxa.combined <- do.call(rbind,taxa.levels) ## concat all levels

taxa.combined <- do.call(rbind,lsTaxa) ## concat all levels
dim(taxa.combined) 

 # 136  72
which(duplicated(rownames(taxa.combined))) 
length(which(duplicated(taxa.combined))) 
# This is due to same feature being collapsed at different levels, 
## so the value of the feature remains the same while the name of taxa changes
## Of the duplicates, keep the one with most characterized levels, so pick from the last of the duplicates
taxa.combined <- taxa.combined[!duplicated(taxa.combined, fromLast = TRUE),] ## genus is last level
length(which(duplicated(taxa.combined))) #0, all dups removed
dim(taxa.combined)
 # 115  72

all.taxa.clr <- as.data.frame(clr_transform(taxa.combined)); dim(all.taxa.clr) 

all.taxa.clr.t <- as.data.frame(t(all.taxa.clr));dim(all.taxa.clr.t)

#save CLR abundances for ileum data
write.table(all.taxa.clr.t, "/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/IBD-eQTL_ileum_16S_Microbial_CLR_abundances.txt", sep="\t", row.names=TRUE, quote=FALSE)

ileum_meta <- sample.info %>% filter(sample.info$biopsy_location=="Ileum" & sample.info$data_type == "biopsy_16S") %>% select(External.ID, Participant.ID)


##order by DNA subject ID
ileum_meta <-ileum_meta[match(rownames(all.taxa.clr.t), ileum_meta$External.ID),]


taxa <- colnames(all.taxa.clr.t)
taxa <- gsub("[", "", taxa, fixed = TRUE)
taxa <- gsub("]", "", taxa, fixed = TRUE)
taxa <- gsub(" ", "_", taxa, fixed = TRUE)
taxa <- gsub(";", ".", taxa, fixed = TRUE)
colnames(all.taxa.clr.t) <- taxa

library(tibble)
all.taxa.clr.t <- rownames_to_column(all.taxa.clr.t, var="External.ID")
rownames(all.taxa.clr.t) <- NULL

taxa_merged <- merge(all.taxa.clr.t, ileum_meta, by = "External.ID")

ileum_taxa <- taxa_merged %>% dplyr::rename(SUBJECT_ID = Participant.ID)

cv_filtered <- cv %>% filter(SUBJECT_ID %in% ileum_taxa$SUBJECT_ID)

cov_merged <- merge(cv_filtered, ileum_taxa, by="SUBJECT_ID")
write.table(cov_merged, file="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/IBD-eQTL_Ileum_covariates_with_microbial_data_10-2-2024.txt",sep="\t",  col.names=T, row.names=F, quote=F)

write.table(taxa, file ="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/1-per-var-list.txt",sep="\t", quote=F, row.names=F,col.names=F)




#make presence/absence covariates

ibd_otu_table <-as.matrix(taxa.combined)


ibd_otu_table <- t(ibd_otu_table)

rows <- rownames(ibd_otu_table)
cols <- colnames(ibd_otu_table)
ibd_otu_table_num <- as.data.frame(matrix(as.numeric(ibd_otu_table), ncol = ncol(ibd_otu_table)))
colnames(ibd_otu_table_num) <- cols
rownames(ibd_otu_table_num) <- rows

##make histograms of abundances to find those with bi-modal distribution
ibd_dichotomous <- ibd_otu_table_num > 0

traits <- c(colnames(ibd_dichotomous))
test <- as.data.frame(ibd_dichotomous) 


min <- nrow(test) / 4
max <-  (nrow(test) / 4) * 3

keep <- colSums(test) > min & colSums(test) < max
ibd_dichotomous <- ibd_dichotomous[,keep]
ibd_dichotomous <- as.data.frame(ibd_dichotomous)
cols <- colnames(ibd_dichotomous)


ibd_dichotomous[,cols] <- lapply(ibd_dichotomous[,cols], as.numeric)
head(ibd_dichotomous)

taxa <- colnames(ibd_dichotomous)
taxa <- gsub("[", "", taxa, fixed = TRUE)
taxa <- gsub("]", "", taxa, fixed = TRUE)
taxa <- gsub(" ", "_", taxa, fixed = TRUE)
taxa <- gsub(";", ".", taxa, fixed = TRUE)
colnames(ibd_dichotomous) <- taxa

ibd_dichotomous <- as.data.frame(ibd_dichotomous)

ibd_dichotomous <- rownames_to_column(ibd_dichotomous, var="External.ID")
rownames(ibd_dichotomous) <- NULL

taxa_merged <- merge(ibd_dichotomous, ileum_meta, by = "External.ID")

ileum_taxa <- taxa_merged %>% dplyr::rename(SUBJECT_ID = Participant.ID)

cv_filtered <- cv %>% filter(SUBJECT_ID %in% ileum_taxa$SUBJECT_ID)

cov_merged <- merge(cv_filtered, ileum_taxa, by="SUBJECT_ID")
write.table(cov_merged, file="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/IBD-eQTL_Ileum_covariates_with_presence_absence_microbial_data_10-2-2024.txt",sep="\t",  col.names=T, row.names=F, quote=F)

write.table(taxa, file ="/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/ileum_microbe_covariates/presence-absence-1-per-var-list.txt",sep="\t", quote=F, row.names=F,col.names=F)


