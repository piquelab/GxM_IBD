##
library(tidyverse)
library(data.table)
library(cowplot)
library(annotables)

outdir <- "./4_INTACT.outs/"

## fn <- "SCAIP_final_bed.gz"
## bed <- fread(fn, header=T, data.table=F)


###############################
### compare results 
###############################

###
### Asthma GBMI_full
traits <- read.table("traits.txt")$V1
trait <- traits[3]
###
fn <- paste(outdir, trait, "_aloft_topPIP_union_intact.txt", sep="")
res_full <- read.table(fn, header=T)

sig <- res_full%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

res2 <- res_full%>%filter(FDR<0.1)



autosome <- as.character(1:22) 
grch38_unq <- grch38%>%
    dplyr::filter(chr%in%autosome, grepl("protein", biotype))%>%
    distinct(ensgene, chr, .keep_all=T)%>%dplyr::select(gene=ensgene, chr, biotype, symbol)

res2 <- res_full%>%left_join(grch38_unq, by=c("Gene"="gene"))

res2%>%filter(symbol%in%c("KIF1B", "TNFRSF14", "NDFIP1"))


###
### old Asthma_UKB
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/enloc_analysis/ALOFT_intact.txt"
old <- read.table(fn, header=T)
sig_old <- old%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

olap <- intersect(sig, sig_old)



#####################################################################
### compare Asthma GBMI_full and the old analysis Asthma UKB   
#####################################################################


traits <- read.table("traits.txt")$V1


###
### GBMI_full
trait <- traits[3]
fn <- paste(outdir, trait, "_aloft_topPIP_union_intact.txt", sep="")
res_full <- read.table(fn, header=T)
res2 <- res_full%>%dplyr::select(Gene, GLCP, zscore=zscore_gwas, pval=pval_gwas)


###
### asthma ukb
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/enloc_analysis/ALOFT_intact.txt"
old <- read.table(fn, header=T)

## twas
fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/Asthma_twas/ALOFT_topPIP_twas.txt.gz"
old_twas <- fread(fn, header=T)%>%dplyr::select(Gene=gene, pval_gwas)

old <- old%>%left_join(old_twas, by="Gene")%>%
   dplyr::select(Gene, GLCP_ukb=GLCP, pval_ukb=pval_gwas) 


### plot data
plotDF <- res2%>%inner_join(old, by="Gene")%>%
    mutate(log10p=-log10(pval), log10p_ukb=-log10(pval_ukb))


###
### compare log10p     
p1 <- ggplot(plotDF, aes(x=log10p_ukb, y=log10p))+
   geom_point(color="lightsteelblue2", size=0.8, shape=1)+
   geom_abline(slope=1, intercept=0, color="grey", linewidth=0.4)+
   xlab(bquote(-log[10]~"("~italic(plain(P))~") from UKB"))+
   ylab(bquote(-log[10]~"("~italic(plain(P))~") from GBMI_full"))+
   ggtitle("ALOFT")+ 
   theme_bw()+
   theme(axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=10))


###
### compare GLCP 
p2 <- ggplot(plotDF, aes(x=GLCP_ukb, y=GLCP))+
   geom_point(color="lightsteelblue2", size=1, shape=1)+
   geom_abline(slope=1, intercept=0, color="grey", linewidth=0.4)+
   xlab("GLCP from UKB")+    
   ylab("GLCP from GBMI_full")+
   ggtitle("ALOFT")+ 
   theme_bw()+
   theme(axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=10))    



####
pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, align="h", axis="tb")
figfn <- paste(outdir, "Figure1_GBMI_fullvsUKB_old.scatter.png", sep="")
ggsave(figfn, pcomb, width=780, height=380, units="px", dpi=120)



################
### overlap
################






###################################
### check new and old
###################################
 
## outdir2 <- "./4_INTACT.outs/check/"

## traits <- read.table("traits.txt")$V1
## trait <- traits[4]
## ###
## fn <- paste(outdir, trait, "_old_aloft_topPIP_union_intact.txt", sep="")
## res <- read.table(fn, header=T)

## sig <- res%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

## ###
## ### old Asthma_UKB
## fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/enloc_analysis/ALOFT_intact.txt"
## old <- read.table(fn, header=T)
## sig_old <- old%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

## olap <- intersect(sig, sig_old)



## ###
## ### plots

## ### UKB_new
## fn <- paste(outdir, trait, "_old_aloft_topPIP_union_intact.txt", sep="")
## res <- read.table(fn, header=T)
## res <- res%>%dplyr::select(Gene, GLCP, zscore=zscore_gwas, pval=pval_gwas)%>%
##     mutate(log10p=-log10(pval))


## ### UKB_old
## fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/enloc_analysis/ALOFT_intact.txt"
## old <- read.table(fn, header=T)

## ## twas
## fn <- "/nfs/rprdata/julong/sc-atac/twas_analysis_2022-10-16/Asthma_twas/ALOFT_topPIP_twas.txt.gz"
## old_twas <- fread(fn, header=T)%>%dplyr::select(Gene=gene, pval_gwas)

## old <- old%>%left_join(old_twas, by="Gene")%>%
##    dplyr::select(Gene, GLCP_ukb=GLCP, pval_ukb=pval_gwas)%>%
##    mutate(log10p_ukb=-log10(pval_ukb))


## ### plot data
## plotDF <- res%>%inner_join(old, by="Gene")


## ###
## ### compare log10p     
## p1 <- ggplot(plotDF, aes(x=log10p_ukb, y=log10p))+
##    geom_point(color="lightsteelblue2", size=0.8, shape=1)+
##    geom_abline(slope=1, intercept=0, color="grey", linewidth=0.4)+
##    xlab(bquote(-log[10]~"("~italic(plain(P))~") from UKB_old"))+
##    ylab(bquote(-log[10]~"("~italic(plain(P))~") from UKB_new"))+
##    ggtitle("ALOFT")+ 
##    theme_bw()+
##    theme(axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))


## ###
## ### compare GLCP 
## p2 <- ggplot(plotDF, aes(x=GLCP_ukb, y=GLCP))+
##    geom_point(color="lightsteelblue2", size=1, shape=1)+
##    geom_abline(slope=1, intercept=0, color="grey", linewidth=0.4)+
##    xlab("GLCP from UKB_old")+    
##    ylab("GLCP from UKB_new")+
##    ggtitle("ALOFT")+ 
##    theme_bw()+
##    theme(axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))    



## ####
## pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, align="h", axis="tb")
## figfn <- paste(outdir2, "Figure1_UKB_NewvsOld.scatter.png", sep="")
## ggsave(figfn, pcomb, width=780, height=380, units="px", dpi=120)
 


## ###
## plotDF <- res2%>%inner_join(old2, by="Gene")

## ### compare GLCP 
## p1 <- ggplot(plotDF, aes(x=GLCP_old, y=GLCP))+
##    geom_point(color="steelblue2", size=0.8, shape=1)+
##    geom_abline(slope=1, intercept=0, color="grey", linewidth=0.4)+
##    xlab("GLCP from old")+
##    ylab("GLCP from current analysis")+
##    ggtitle(paste(trait0, "(ALOFT)"))+ 
##    theme_bw()+
##    theme(axis.text=element_text(size=10),
##          axis.title=element_text(size=10),
##          plot.title=element_text(hjust=0.5, size=10))    

## ###
## figfn <- paste(outdir2, "Figure0.1_", trait, "_new_old.GLCP.scatter.png", sep="")
## ggsave(figfn, p1, width=380, height=380, units="px", dpi=120)
     
## ## p2 <- ggplot(plotDF, aes(x=zscore_old, y=zscore_gwas))+
## ##    geom_point(color="steelblue2", size=0.8, shape=1)+
## ##    geom_abline(slope=1, intercept=0, color="grey", linewidth=0.4)+
## ##    xlab("Zscore from old")+
## ##    ylab("Zscore from current analysis")+    
## ##    ggtitle(paste(trait0, "(ALOFT)"))+ 
## ##    theme_bw()+
## ##    theme(axis.text=element_text(size=10),
## ##          axis.title=element_text(size=10),
## ##          plot.title=element_text(hjust=0.5, size=10))

## ####
## pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, align="h", axis="tb")
## figfn <- paste(outdir2, "Figure0_", trait, "_new_old.scatter.png", sep="")
## ggsave(figfn, pcomb, width=780, height=380, units="px", dpi=120)

 

## ###
## ###

## sig_new <- res%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

## sig_old <- old%>%filter(FDR<0.1)%>%pull(Gene)%>%unique()

## length(intersect(sig_new, sig_old))
