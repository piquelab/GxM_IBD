#
library(tidyverse)
library(data.table)
library(annotables)

options(scipen=16)


rm(list=ls())



outdir <- "./torus_input/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


####################
### torus format ###
####################

tissue <- "IBD_ileum"
infn_eqtl <- "../0_FastQTL/fastQTL.output/IBD_ileum_nominals_all_chunks.txt.gz"
opfn <- gzfile(paste(outdir, tissue, ".eQTL.txt.gz", sep=""))

###
### eQTL summary results
res <- fread(infn_eqtl, header=F, data.table=F, stringsAsFactors=F)

zscore <- abs(qnorm(res$V4/2))*sign(res$V5)

summ <- data.frame(SNP=res$V2, gene=res$V1, beta=res$V5, "t-stat"=zscore, "p-value"=res$V4)
  
###output 
write.table(summ, opfn, quote=F, row.names=F)
close(opfn)


###
geneList <- unique(res$V1)
opfn_gene <- "geneList.txt"
write.table(geneList, file=opfn_gene, quote=F, row.names=F, col.names=F)

##
snpList <- unique(res$V2)
opfn_snp <- "snpList.txt"
write.table(snpList, file=opfn_snp, quote=F, row.names=F, col.names=F)




##########################
### gene map for torus ###
##########################

###
### annotation files
## anno <- read.table("/wsu/home/groups/piquelab/data/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gff3.gz", header=F, stringsAsFactors=F) ##grch37

###
### gene annotation files from grch38
chr_auto <- as.character(1:22)
anno <- grch38%>%dplyr::filter(chr%in%chr_auto, grepl("protein_coding", biotype))%>%
    dplyr::select(ensgene, chr, s0=start, s1=end, strand)

anno2 <- anno%>%
   mutate(start=ifelse(strand=="+", s0, s1),
          end=ifelse(strand=="-", s1, s0))%>%
   dplyr::select(-s0, -s1, -strand) 


###
### gene used for Rectum 
geneList <- read.table("geneList.txt")$V1 ### 20,315 gene

anno3 <- anno2%>%dplyr::filter(ensgene%in%geneList)%>%
   dplyr::select(gene=ensgene, chr, start) 
anno3$start2 <- anno3$start


###
opfn <- gzfile(paste(outdir, "zzz_gene.map.gz", sep=""))
write.table(anno3, opfn, sep="\t", quote=F, row.names=F, col.names=F)

 


#########################################
### SNP map file and annotation files ###
#########################################

fn <- "../snpinfor.txt.gz"
DF_snp <- fread(fn, header=F, data.table=F)

DF2 <- DF_snp%>%dplyr::select(SNP=V5, chr=V1, pos=V2)

### output all snp map
opfn2 <- gzfile(paste(outdir, "zzz_snp.map.gz", sep=""))
write.table(DF2, opfn2,  quote=F, row.names=F, col.names=F)



