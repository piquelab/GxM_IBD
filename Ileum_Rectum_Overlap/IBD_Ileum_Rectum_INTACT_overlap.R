library(tidyverse)
library(ggplot2)
library(data.table)
library(stringr)
library(annotables)
library(RColorBrewer)
library(Hmisc)

#overlap of rectum and ileum INTACT results

dataDir <- "/rs/rs_grp_ibdeqtl/IBD_shreya_2024/IBD_ileum_2024-10-22/4_TWAS_smr/INTACT_output/eqtl_topPIP"

traits <- c("IBD.EUR.Crohns_Disease", "IBD.EUR.Inflammatory_Bowel_Disease", "IBD.EUR.Ulcerative_Colitis")

resList <- list()
#read in files and extract signficiant risk associated eGenes

for(i in 1:length(traits)) {
	trait <- traits[i]
	fname <- paste0(dataDir, "/", trait, "_IBD_ileum_intact.txt")
	res <- read.table(fname, stringsAsFactors=F, header=T, comment="")

			res <- res %>% select(Gene, GRCP, GLCP, zscore, LFDR, FDR)
			res$Trait <- trait


			resList[[i]] <- res

		}


intact_all <- data.frame(matrix(ncol = ncol(resList[[1]]), nrow = 0))
        colnames(intact_all) <- colnames(resList[[1]])

        for(n in 1:length(resList)) {
            intact_all <- rbind(intact_all, resList[[n]])
        }


intact_sig <- intact_all %>% filter(FDR < 0.1)
length(intact_sig$Gene) 
#209
length(unique(intact_sig$Gene))

ileum_intact_sig <- intact_sig

dataDir <- "/rs/rs_grp_ibdeqtl/IBD_shreya_2024/IBD_Rectum2_2024-04-14/4_TWAS_smr/INTACT_output/eqtl_topPIP"

traits <- c("IBD.EUR.Crohns_Disease", "IBD.EUR.Inflammatory_Bowel_Disease", "IBD.EUR.Ulcerative_Colitis")

resList <- list()
#read in files and extract signficiant risk associated eGenes

for(i in 1:length(traits)) {
	trait <- traits[i]
	fname <- paste0(dataDir, "/", trait, "_IBD_Rectum_intact.txt")
	res <- read.table(fname, stringsAsFactors=F, header=T, comment="")

			res <- res %>% select(Gene, GRCP, GLCP, zscore, LFDR, FDR)
			res$Trait <- trait


			resList[[i]] <- res

		}


rectum_intact_all <- data.frame(matrix(ncol = ncol(resList[[1]]), nrow = 0))
        colnames(rectum_intact_all) <- colnames(resList[[1]])

        for(n in 1:length(resList)) {
            rectum_intact_all <- rbind(rectum_intact_all, resList[[n]])
        }


rectum_intact_sig <- rectum_intact_all %>% filter(FDR < 0.1)

rectum_risk_genes <- unique(rectum_intact_sig$Gene)

ileum_risk_genes <- unique(ileum_intact_sig$Gene)

table(ileum_risk_genes %in% rectum_risk_genes)


# FALSE  TRUE
#    62    69
table(rectum_risk_genes %in% ileum_risk_genes)

 table(rectum_risk_genes %in% ileum_risk_genes)

# FALSE  TRUE
#    53    69

#174 unique genes across both body sites

#69 shared

#make barplot showing overlap of ileum and rectum INTACT genes

Body_sites <- c(rep("Rectum", 3), rep("Ileum", 3))
significance <- rep(c("IBD Risk Gene Rectum", "Shared IBD Risk Gene", "IBD Risk Gene Ileum"), 2)
numgenes <- c(53, 69, 0, 0, 69, 62)

df <- data.frame(Body_sites, significance, numgenes)
df$significance <- factor(df$significance, levels = c("IBD Risk Gene Rectum", "Shared IBD Risk Gene", "IBD Risk Gene Ileum"))

png("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD-eQTL_BodySite_INTACT_overlap_barplot.png")
p <- ggplot(df, aes(fill = significance, y = numgenes, x = Body_sites)) +
		geom_bar(position="stack", stat="identity", width = 0.5) +
     	theme_classic() +
      	scale_fill_brewer(palette = "Dark2") +
      	ylab("# of IBD Risk genes (INTACT FDR = 10%)") +
      	coord_flip() +
      	theme(axis.title.y = element_blank())
print(p)
dev.off()

png("/rs/rs_grp_ibdeqtl/FastQTL/ileum_protein_coding_residuals/qnorm_residual_analysis/GxM_INTACT/IBD-eQTL_BodySite_INTACT_overlap_barplot_no_legend.png", width = 600, height = 200)
p <- ggplot(df, aes(fill = significance, y = numgenes, x = Body_sites)) +
		geom_bar(position="stack", stat="identity", width = 0.5) +
     	theme_classic() +
      	scale_fill_brewer(palette = "Dark2") +
      	ylab("# of IBD Risk genes (INTACT FDR = 10%)") +
      	coord_flip() +
      	theme(legend.position = "none", axis.title.y = element_blank())
print(p)
dev.off()

