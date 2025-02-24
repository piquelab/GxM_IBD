###
###
library(tidyverse)
library(data.table)
library(cowplot)
library(annotables)

outdir <- "./INTACT_output/check/"

if ( ! file.exists(outdir)) dir.create(outdir, recursive=T)


########################
### testing genes 
#########################

###
### new data

##
fn <- "./eQTL_results/IBD_Rectum_eqtl.txt"
fast <- read.table(fn, header=T)
##
fn <- "./eQTL_results/IBD_Rectum_PIP.txt"
dap <- read.table(fn, header=T)

gene <- union(fast$gene, dap$gene)


####
### old data
fn <- "/rs/rs_grp_ibdeqtl/IBD_shreya_2024/IBD_Rectum/3_gwas_prepare/eQTL_results/IBD_Rectum_eqtl.txt.gz"
fast_old <- fread(fn, header=T, data.table=F)

fn <- "/rs/rs_grp_ibdeqtl/IBD_shreya_2024/IBD_Rectum/3_gwas_prepare/eQTL_results/IBD_Rectum_PIP.txt.gz"
dap_old <- fread(fn, header=T, data.table=F)

gene_old <- union(fast_old$gene, dap_old$gene)



#############################
### compare new and old 
#############################



twas_appro <- c("minP", "topPIP")
traits <- read.table("traits_ls.txt")$V1




for (mm in twas_appro){
for (trait in traits){
    
###
fn <- paste("./INTACT_output/eqtl_", mm, "/",  trait, "_IBD_Rectum_intact.txt", sep="")
res <- read.table(fn, header=T)
res2 <- res%>%dplyr::select(Gene, GLCP, zscore, pval=pval_gwas)


###
### old
    prefix <- "/rs/rs_grp_ibdeqtl/IBD_shreya_2024/IBD_Rectum/4_TWAS_smr/"    
    fn2 <- paste(prefix, fn, sep="") 
    old <- read.table(fn2, header=T)

    fn0 <- paste(prefix, "1_SMR_output/", trait, "_", mm, "_twas.txt.gz", sep="")
    twas <- fread(fn0, header=T, data.table=F)%>%dplyr::select(Gene=gene, pval_gwas)
    old <- old%>%left_join(twas, by="Gene")
    old <- old%>%dplyr::select(Gene, GLCP_old=GLCP, pval_old=pval_gwas)
    
###
### plot data
    plotDF <- res2%>%inner_join(old, by="Gene")%>%
       mutate(log10p=-log10(pval), log10p_old=-log10(pval_old))
   comb <- paste(mm, gsub("^IBD.EUR.", "", trait))
    
###
### compare log10p     
p1 <- ggplot(plotDF, aes(x=log10p_old, y=log10p))+
   geom_point(color="lightsteelblue2", size=0.8, shape=1)+
   geom_abline(slope=1, intercept=0, color="grey", linewidth=0.4)+
   xlab(bquote(-log[10]~"("~italic(plain(P))~") from old"))+
   ylab(bquote(-log[10]~"("~italic(plain(P))~") from new"))+
   ggtitle(comb)+ 
   theme_bw()+
   theme(axis.text=element_text(size=10),
         axis.title=element_text(size=10),
         plot.title=element_text(hjust=0.5, size=10))


 ###
 ### compare GLCP 
 p2 <- ggplot(plotDF, aes(x=GLCP_old, y=GLCP))+
    geom_point(color="lightsteelblue2", size=1, shape=1)+
    geom_abline(slope=1, intercept=0, color="grey", linewidth=0.4)+
    xlab("GLCP from old")+    
    ylab("GLCP from new")+
    ggtitle(comb)+ 
    theme_bw()+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title=element_text(hjust=0.5, size=10))    

cat(comb, "\n")

####
pcomb <- plot_grid(p1, p2, nrow=1, ncol=2, align="h", axis="tb")
figfn <- paste(outdir, "Figure1_",  mm, "_", trait, "_old.scatter.png", sep="")
ggsave(figfn, pcomb, width=780, height=380, units="px", dpi=120)

}
} 
    

    
