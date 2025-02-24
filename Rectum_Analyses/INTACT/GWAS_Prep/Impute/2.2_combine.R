###
###
library(tidyverse)
library(data.table)
options(scipen=15)

rm(list=ls())


###
### directory for outputs
outdir <- "./2_impute.outs/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)
 


traits <- read.table("traits_ls.txt")$V1

###
### loop for trait
for ( trait in traits){
    
###
### gwas data
fn <- paste("../gwas_data/", trait, "_gwas.txt.gz", sep="")
summ <- fread(fn, header=T, data.table=F)
summ <- summ%>%mutate(id_b38_0=id_b38)
    
### rbind

## summ2 <- rbind(summ, tmp)
## ###
## gfn <- paste(outdir, trait, "_gwas_imput.txt.gz", sep="")
## fwrite(summ2, file=gfn, sep=" ", quote=F, row.names=F, na=NA)

## ## ### missing
## gfn2 <- paste(outdir, trait, "_gwas_miss.txt.gz", sep="")
## fwrite(tmp, gfn2, sep=" ", quote=F, row.names=F, na=NA)

fn2 <- paste(outdir, trait, "/all_gwas_impute.txt.gz", sep="")
x2 <- fread(fn2, header=F, data.table=F)
names(x2) <- c("id_b38", "chr", "pos", "zscore", "pval", "id_b38_0")
    
summ2 <- rbind(summ, x2)

gfn <- paste(outdir, trait, "_gwas_imput.txt.gz", sep="")
fwrite(summ2, file=gfn, sep=" ", quote=F, row.names=F, na=NA)

cat(trait, "\n")

}    
    
###
### END
