for j in $(seq 1 30); do head -n $(($j+1)) IBD_eQTL_Ileum-PCcovariates-FastQTL.txt > covariates/IBD_eQTL_Ileum-PC1-$j.full.covariates-FastQTL.txt; done


head -n 1 IBD_eQTL_Ileum-PCcovariates-FastQTL.txt > covariates/IBD_eQTL_Ileum-PC1-0.full.covariates-FastQTL.txt