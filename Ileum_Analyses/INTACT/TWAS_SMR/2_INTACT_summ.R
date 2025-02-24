##
library(colorspace) ##, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library")
library(tidyverse) ###, lib.loc="/wsu/el7/groups/piquelab/R/4.1.0/lib64/R/library")
library(data.table)
library(annotables)
library(INTACT, lib.loc="/wsu/home/ha/ha21/ha2164/Bin/Rpackages/")
library(openxlsx)

rm(list=ls())




#####################################################
### run INTACT to combine colocalzaition and TWAS ###
#####################################################


traits <- sort(read.table("traits_ls.txt")$V1) 
twas_appro <- c("minP", "topPIP")
tissue <- "IBD_ileum"

###
for ( ii in twas_appro){

##    
outdir2 <- paste("./INTACT_output/eqtl_", ii, "/", sep="")
if ( !file.exists(outdir2)) dir.create(outdir2, recursive=T)

for ( trait in traits){

cat(ii, trait, tissue, ",")
    
##    
fn <- paste("../3_enloc/enloc_output/", trait, "_", tissue, ".enloc.gene.out", sep="")
enloc <- read.table(fn, header=T)

    
###
### smr
fn <- paste("./1_SMR_output/", trait, "_", ii, "_twas.txt.gz", sep="")
twas <- fread(fn, header=T, data.table=F)
names(twas)[1] <- "Gene"

## combine enloc and z-score of twas
DF <- enloc%>%inner_join(twas, by="Gene")%>%drop_na(zscore)
    
## intact analysis
res_intact <- intact(GLCP=DF$GLCP, z_vec=DF$zscore)
DF$PCG <- res_intact
### FDR 
DF2 <- DF%>%mutate(LFDR=1-PCG)%>%arrange(LFDR)
x <- DF2$LFDR
FDR <- cumsum(x)/1:length(x)
DF2$FDR <- FDR

### 
opfn <- paste(outdir2, trait, "_", tissue,  "_intact.txt", sep="")
write.table(DF2, opfn, row.names=F, quote=F, col.names=T)

cat("SMR eqtl:", sum(FDR<0.05), "\n")    

} ### END loop for traits
} ### END loop for TWAS approach    







#######################
### summary results ###
#######################

## rm(list=ls())

## traits <- sort(read.table("traits_ls.txt")$V1) 
## tissues <- c("Colon_Sigmoid", "Colon_Transverse")
## twas <- c("predixcan", "ptwas")

## outdir <- "./INTACT_output/"

## summ <- NULL
## for (trait in traits){
##    ###
##    for (tissue in tissues){
##        ###
##        nums <- map_dbl(twas, function(ii){
##            ##
##            fn <- paste(outdir, trait, "/", trait, "_", tissue, "_", ii, "_intact.txt", sep="")
##            res <- read.table(fn, header=T)%>%filter(FDR<0.05)
##            nrow(res)
##        })
##       df2 <- data.frame(traits=rep(trait, 2), tissues=rep(tissue, 2), twas_approach=twas, nsig_0.05=nums)
##       summ <- rbind(summ, df2)
##    }
## }

## opfn <- paste(outdir, "IBD_summary.xlsx", sep="")
## write.xlsx(summ, file=opfn, overwrite=T)
    
       
           
