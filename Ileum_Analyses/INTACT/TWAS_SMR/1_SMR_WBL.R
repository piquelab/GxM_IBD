##
library(tidyverse)
library(data.table)

rm(list=ls())
 
###
outdir <- "./1_SMR_output/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


### passing argument
args=commandArgs(trailingOnly=T)
if (length(args)>0){
   trait <- args[1]
 }else{
   trait <- "IBD.EUR.Crohns_Disease"
}

tissue <- "IBD_ileum"
## traits <- read.table("traits_of_interest.txt")$V1
## traits <- sort(traits)
 

###
### gwas data _gwas_imput.txt.gz
fn <- paste0("./gwas_imputefile/", trait, "_gwas_imput.txt.gz")
summ <- fread(fn, header=T, data.table=F)
summ2 <- summ%>%distinct(id_b38_0, .keep_all=T)
p_gwas <- summ2$pval
names(p_gwas) <- summ2$id_b38_0
zval <- summ2$zscore
names(zval) <- summ2$id_b38_0


###
### option-1, fastqtl, use minimum pvalue to pick gwas pvalue

fn <- paste("./eQTL_results/", tissue, "_eqtl.txt", sep="")
fast <- fread(fn, header=T, data.table=F)
fast <- fast%>%mutate(pval_gwas=p_gwas[id_b38], zscore=zval[id_b38])
 
res_gene <- fast%>%drop_na(pval_gwas)%>%
    group_by(gene)%>%slice_min(order_by=pval, n=1)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%
    ungroup()%>%as.data.frame()

###
gfn <- gzfile(paste(outdir, trait, "_minP_twas.txt.gz", sep=""))
write.table(res_gene, gfn, row.names=F, col.names=T, quote=F)




###
### option-2, dap-g, use top PIP to pick gwas pvalue 
fn <- paste("./eQTL_results/", tissue, "_PIP.txt", sep="")
dap <- fread(fn, header=T, data.table=F)
dap <- dap%>%mutate(pval_gwas=p_gwas[id_b38], zscore=zval[id_b38])

res2_gene <- dap%>%drop_na(pval_gwas)%>%
    group_by(gene)%>%slice_max(order_by=PIP, n=1)%>%
    slice_min(order_by=pval_gwas, n=1, with_ties=F)%>%
    ungroup()%>%as.data.frame()

###
gfn2 <- gzfile(paste(outdir, trait, "_topPIP_twas.txt.gz", sep=""))
write.table(res2_gene, gfn2, row.names=F, col.names=T, quote=F)

###
###END
    
    
