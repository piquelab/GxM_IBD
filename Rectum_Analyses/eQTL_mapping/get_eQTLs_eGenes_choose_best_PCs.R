## this script read in all 0-20 PCs permutations runs to get the # of PCs maximizing eGenes:
# 8/26/2019 JR

##Adpated by Shreya Nirmalan for IBD eQTL project 2-22-23 - used for rectum protein coding residuals 11-19-23
library(data.table)
library(qqman)
library(qvalue)

FDR <- 0.1

## load all 21 permutation files:
# make a list of all files and read them in:
## files.list <- list.files(path="../permutations/results/", pattern = "*.permutations.eQTL.txt.gz")
all_PCs <- lapply(Sys.glob("../output/*.permutations.eQTL.txt.gz"), fread)
names(all_PCs) <- gsub(".*-(.*).permutations.eQTL.txt.gz", "\\1", Sys.glob("../output/*.permutations.eQTL.txt.gz"))
# add colnames:
for(i in 1:length(all_PCs)){
    colnames(all_PCs[[i]]) <- c("pid", "nvar", "shape1", "shape2", "dummy", "sid", "dist", "npval", "slope", "ppval", "bpval")
}
# add qvalue to each:
res <- data.frame("PCs"=integer(0),"eGenes"=integer(0))
for(i in 1:length(all_PCs)){
all_PCs[[i]]$bqval <- qvalue(all_PCs[[i]]$bpval)$qvalues
res <- rbind(res, as.numeric(c(names(all_PCs)[i],sum(all_PCs[[i]]$bqval<FDR,na.rm =TRUE)))
)    }
colnames(res) <- c("PCs","eGenes")
best.index <- which(res$eGenes==max(res$eGenes))
res <- res[order(res$PCs),]
 best.PCs <- res[res$eGenes==max(res$eGenes),"PCs"]

# chose best # of PCs:
res[which(res$eGenes==max(res$eGenes)),]
#PC 19 3455 eGenes

# save:
write.table(res, "eGenes-per-GEPCs.txt", sep='\t', quote=F, row.names=F)

# save the best results:
write.table(all_PCs[[best.index]], paste0("../results/FastQTL_results_best_", best.PCs, "GEPCs.txt"), sep='\t', quote=F, row.names=F)

# subset to significant only:
pc19_signif_pairs <- all_PCs[[best.index]]
pc19_signif_pairs <- pc19_signif_pairs[pc19_signif_pairs$bqval<FDR,]

pairs <- pc19_signif_pairs[,c(1,6),]
write.table(pairs, file="./PC19_significant_topeeQTL_pairs.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

##save SNP IDs
snps <- pc19_signif_pairs[,"sid"]
write.table(snps, file="./PC19_significant_topeeQTL_snps.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# # split the variant info into multiple columns:
# library(tidyr)
# pc20_signif_pairs <- separate(pc20_signif_pairs,6,into=c(NA, "rsID"),sep=";", remove=FALSE)
# pc20_signif_pairs <- separate(pc20_signif_pairs,6,into=c("chr","position","var1","var2"),sep=":", remove=FALSE)
# write.table(pc20_signif_pairs, file="./PC1-20.permutations.eQTL_signif.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# # save just the chr, pos:
# write.table(pc20_signif_pairs[,c("chr","position")], file="./PC20_significant_topeeQTL_positions.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
# # save chr, SNP pairs:
# pairs <- pc20_signif_pairs[,c(1,6),]
# write.table(pairs, file="./PC20_significant_topeeQTL_pairs.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# ### END 10/19/2020


##make qqplot with both the nominal p-value and permutation p-value for PC20
pc19_results <- all_PCs[[best.index]]

library(ggplot2)
library(tidyverse)
ci=0.95


pc19_results_bp <- pc19_results %>% select(pid, sid, bpval) %>% filter(!is.na(bpval)) %>%
            arrange(bpval) %>%
            mutate(r=rank(bpval, ties.method = "random"),
                   pexp=r/length(bpval),
                   clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = r, shape2 = length(bpval)-r)),
                   cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = r, shape2 = length(bpval)-r)))

pc19_results_bp$group <- "Corrected p-value"
pc19_results_bp <- pc19_results_bp %>% dplyr::rename(pval = bpval) 



                pdf("IBD-eQTL_Rectum_PC19_eGene_qqplotpermuted_pvalue_only_no_legend.pdf")          
        p1 <- ggplot(pc19_results_bp, aes(x=-log10(pexp),y=-log10(pval))) +
            geom_ribbon(mapping = aes(x = -log10(pexp), ymin = clower, ymax = cupper),
              alpha = 0.1,color="darkgray") +
            geom_point() +
            geom_abline(slope=1,intercept=0) +
        ##    facet_grid(Origin ~ Location) +
            xlab(expression(Expected -log[10](p))) +
            ylab(expression(Observed -log[10](p))) + 
            ggtitle(paste0("IBD eQTL Rectum PC19 eGene QQ Plot")) +
            theme_classic() +
            theme(legend.title= element_blank(), axis.title.x = element_text(size = rel(1.2)), axis.title.y = element_text(size = rel(1.2)), legend.text = element_blank(), plot.title = element_text(hjust=0.5,size = rel(1.3)))
        print(p1)
        dev.off()