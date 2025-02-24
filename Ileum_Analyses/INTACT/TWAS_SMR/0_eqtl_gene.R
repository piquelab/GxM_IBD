###
library(tidyverse)
library(data.table)

##
options(scipen=16)

rm(list=ls())

##
outdir <- "./eQTL_results/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=F)

###
###
### SNP vcf
fn <- "../snpinfor.txt.gz"
DF_snp <- fread(fn, header=F, data.table=F)%>%filter(V1%in%as.character(1:22))
DF_snp2 <- DF_snp%>%mutate(id_b38=paste0("chr", paste(V1, V2, V3, V4, "b38", sep="_")))%>%
    dplyr::select(id_b38, rs=V5)

DF_snp2 <- DF_snp2%>%distinct(rs, .keep_all=T)

id_chr_pos <- DF_snp2$id_b38
names(id_chr_pos) <- DF_snp2$rs

####
#### SNPs in eqtls
tissue <- "IBD_ileum"
infn_eqtl <- paste("../0_FastQTL/fastQTL.output/", tissue, "_nominals_all_chunks.txt.gz", sep="")
opfn_eqtl <- paste(outdir, tissue, "_eqtl.txt", sep="")
n0 <- 10


### 
res <- fread(infn_eqtl, header=F, data.table=F, stringsAsFactors=F)
names(res) <- c("gene", "rs", "DTSS", "pval", "beta")
 
res2 <- res%>%group_by(gene)%>%slice_min(order_by=pval, n=n0, with_ties=T)%>%as.data.frame()
 
res3 <- res2%>%filter(rs%in%DF_snp2$rs)%>%mutate(id_b38=id_chr_pos[rs])

### outputs
write.table(res3, file=opfn_eqtl, row.names=F, col.names=T, quote=F)




###
### SNPs in PIP
tissue <- "IBD_ileum"
infn_dap <- paste("../1_DAP-G/", tissue, "_PIP.txt.gz", sep="")
opfn_dap <- paste(outdir, tissue, "_PIP.txt", sep="")
n0 <- 10

dap <- fread(infn_dap, header=F, data.table=F, stringsAsFactors=F)
names(dap) <- c("gene", "rs", "PIP")

dap2 <- dap%>%group_by(gene)%>%slice_max(order_by=PIP, n=n0, with_ties=T)%>%as.data.frame()

dap3 <- dap2%>%filter(rs%in%DF_snp2$rs)%>%mutate(id_b38=id_chr_pos[rs])

### outputs

write.table(dap3, file=opfn_dap, row.names=F, col.names=T, quote=F)


###
### END

