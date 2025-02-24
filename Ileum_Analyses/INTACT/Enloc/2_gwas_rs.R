###
###
library(tidyverse)
library(data.table)

##
outdir <- "./gwas_PIP2/"
if ( !file.exists(outdir) ) dir.create(outdir, showWarnings=F, recursive=T)


### replace SNP id and keep the same to those in eqtl mapping.




df_snp <- fread("../snpinfor.txt.gz", header=F, data.table=F)
df2 <- df_snp%>%mutate(id_b38=paste0("chr", paste(V1, V2, V3, V4, "b38", sep="_")))
rs_id <- df2$V5
names(rs_id) <- df2$id_b38


###
traits <- read.table("./traits_ls.txt")$V1
for (ii in traits){

###
cat(ii, "\n")
    
fn <- paste("./gwas_PIP/", ii, ".pip.gz", sep="")
x <- fread(fn, header=F, data.table=F)    

x2 <- x%>%mutate(is_rs=V1%in%df2$id_b38, rs=ifelse(is_rs, rs_id[V1], V1))

### output    
x3 <- x2%>%dplyr::select(rs, V2, V3, V4)
opfn <- paste0(outdir, ii, ".pip.gz")
fwrite(x3, file=opfn, row.names=F, col.names=F, sep=" ")
}

  
