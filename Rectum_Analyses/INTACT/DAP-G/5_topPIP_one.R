##
library(Matrix)
library(tidyverse)
library(data.table)
## library(clusterProfiler)
## library(org.Hs.eg.db)
## library(annotables)
## library(ggplot2)
## library(cowplot)
## library(RColorBrewer)

rm(list=ls())


###parsing arguments
args=commandArgs(trailingOnly=T)
if ( length(args)>0){
   ###
   geneFile <- args[1]
}else{
    geneFile <- "splitGene000"
}



###
outdir <- "./5_summary.outs/IBD_Rectum/"
if ( !file.exists(outdir)) dir.create(outdir, showWarnings=F, recursive=T)


## peak SNP with maximum PIP 

###
###
fn <- paste("./geneList/", geneFile, sep="")
geneList <- read.table(fn)$V1
dap <- NULL
for (ens in geneList){
   ###
   cat(ens, "\n") 
   dapfn <- paste("./dap-g_parse_outs/", ens, ".SNP.out", sep="")
   if ( file.exists(dapfn)&file.size(dapfn)>0){
      ###    
      res2 <- read.table(dapfn)%>%mutate(gene=ens)%>%dplyr::select(gene, SNP=V2, PIP=V3)
      dap <- rbind(dap, res2)
  }
}    
 
opfn <- paste(outdir, "dap_", geneFile, ".txt", sep="")
write.table(dap, file=opfn, sep="\t", row.names=F, quote=F, col.names=F)

###
