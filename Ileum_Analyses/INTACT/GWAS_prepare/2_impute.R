###
###
library(tidyverse)
library(data.table)
options(scipen=15)

rm(list=ls())

### parsing argument
args=commandArgs(trailingOnly=T)
if ( length(args)>0){
   ##
   trait <- args[1] 
   snpfn <- args[2]
}else{
   trait <- "IBD.EUR.Crohns_Disease"
   snpfn <- "zzz_splitSNP0000"
}


###
### directory for outputs
outdir <-paste("./2_impute.outs/", trait, "/", sep="")
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
 


###
### gwas data
fn <- paste("../gwas_data/", trait, "_gwas.txt.gz", sep="")
summ <- fread(fn, header=T, data.table=F)



###
### missing SNP files
fn <- paste("./1_missing.outs/", trait, "/", snpfn, sep="")
snp <- read.table(fn, header=F)$V1


##
x <- str_split(snp, "_", simplify=T)
bed <- data.frame(id_b38=snp, chr=x[,1], pos=x[,2])

### missing SNPs, 3499


###
### imputation
nsnp <- nrow(bed)
time0 <- Sys.time()
tmp <- lapply(1:nsnp, function(i){
   ##
   if (i%%10==0) cat(i,"\n") 
   chr_i <- bed$chr[i]
   pos_i <- as.integer(bed$pos[i])

   ## +/- 10 kb 
   s0 <- pos_i-1e+04
   s1 <- pos_i+1e+04
    
   tmp2 <- summ%>%filter(chr==chr_i)
   pos <- tmp2$pos
   sub0 <- pos>=s0&pos<s1 
   nsnp <- sum(sub0)
    
   if ( nsnp>0){
      tmp2 <- tmp2[sub0,] 
      dtss <- abs(as.numeric(tmp2$pos)-pos_i)
      imin <- which.min(dtss)
      tmp3 <- tmp2[imin,]
      tmp3$id_b38_0 <- bed$id_b38[i]
   }else{
      tmp3 <- NULL
   }   
   tmp3
})
tmp <- do.call(rbind, tmp)


##############
### output ###
##############

### rbind
## summ <- summ%>%mutate(id_b38_0=id_b38)
## summ2 <- rbind(summ, tmp)
## ###
## gfn <- paste(outdir, trait, "_gwas_imput.txt.gz", sep="")
## fwrite(summ2, file=gfn, sep=" ", quote=F, row.names=F, na=NA)

## ## ### missing
## gfn2 <- paste(outdir, trait, "_gwas_miss.txt.gz", sep="")
## fwrite(tmp, gfn2, sep=" ", quote=F, row.names=F, na=NA)

gfn2 <- paste(outdir, snpfn, "_miss.txt.gz", sep="")
fwrite(tmp, gfn2, sep=" ", quote=F, row.names=F, na=NA)

time1 <- Sys.time()
diff0 <- difftime(time1, time0, units="secs")
cat(diff0, "\n")


### END
