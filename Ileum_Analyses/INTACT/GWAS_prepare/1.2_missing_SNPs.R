
##
library(tidyverse)
library(data.table)

options(scipen=16)
 
rm(list=ls())


outdir <- "./1_missing.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)



###
### SNP vcf
prefix <- "/rs/rs_grp_ibdeqtl/IBD_shreya_2024/"

dir_ls <- c("./IBD_Rectum2_2024-04-14/", "./IBD_ileum_2024-10-22/")


DF_snp <- map_dfr(dir_ls, function(ii){
    ##
    fn0 <- paste(prefix, ii, "snpinfor.txt.gz", sep="")            
    DF_snp <- fread(fn0, header=F, data.table=F)%>%filter(V1%in%as.character(1:22))
    DF_snp2 <- DF_snp%>%mutate(id_b38=paste0("chr", paste(V1, V2, V3, V4, "b38", sep="_")))%>%
       dplyr::select(id_b38, rs=V5)
    ##
    DF_snp2 <- DF_snp2%>%distinct(rs, .keep_all=T)
    DF_snp2
})

DF_snp2 <- DF_snp%>%distinct(rs, .keep_all=T)
 

####
#### SNPs in eqtls
tissue_ls <- c("IBD_Rectum", "IBD_ileum")

snpSel <- lapply(tissue_ls, function(ii){
   ###
   fn0 <- paste(ii, ".snpList.txt", sep="")
   snp0 <- read.table(fn0, header=F)$V1
   snp0
})%>%unlist()%>%unique()    




DF_snp2 <- DF_snp2%>%filter(rs%in%snpSel)




###
###
traits <- read.table("traits_ls.txt")$V1
traits <- sort(traits)

for (ii in traits){
   ### 
   fn <- paste("../gwas_data/", ii, "_gwas.txt.gz", sep="")
   summ <- fread(fn, header=T, data.table=F)


   snp_miss <- DF_snp2%>%filter(!id_b38%in%summ$id_b38)

   cat(ii, nrow(snp_miss), "\n") 
   ### output 
   opfn <- paste(outdir, ii, "_missing.txt", sep="")
   write.table(snp_miss, opfn, quote=F, row.names=F, col.names=T)
}





##################
### imputation ###
##################

## rm(list=ls())

## outdir2 <- "./gwas_imputefile/"


## fn <- paste(outdir2, "snp_missing.txt", sep="")
## snp_miss <- read.table(fn)$V1


## traits <- read.table("traits_ls.txt")$V1
## traits <- sort(traits)

## for (ii in traits){
##    ##
##    ii <- traits[1]
##    fn <- paste("./gwas_data/", ii, "_gwas.txt.gz", sep="")
##    x <- fread(fn, header=T, data.table=F)

##    snp2 <- snp_miss[!snp_miss%in%x$id_b38]
   
##    chr_i <- bed$chr[i]
##    pos_i <- as.integer(bed$pos[i])

##    ## +/- 10 kb 
##    s0 <- pos_i-1e+04
##    s1 <- pos_i+1e+04
    
##    tmp2 <- summ%>%filter(chr==chr_i)
##    pos <- tmp2$pos
##    sub0 <- pos>=s0&pos<s1 
##    nsnp <- sum(sub0)
    
##    if ( nsnp>0){
##       tmp2 <- tmp2[sub0,] 
##       dtss <- abs(as.numeric(tmp2$pos)-pos_i)
##       imin <- which.min(dtss)
##       tmp3 <- data.frame(id_b38=bed$id_b38[i], chr=chr_i, pos=pos_i,
##                          zscore=tmp2$zscore[imin], pval=tmp2$pval[imin],
##                          id_b38_close=tmp2$id_b38[imin])
##    }else{
##       tmp3 <- NULL
##    }   
##    tmp3
## })


## time1 <- Sys.time()
## diff0 <- difftime(time1, time0, units="secs")
## cat(diff0, "\n")    
