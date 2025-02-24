#!/bin/bash/

cd $PWD

cat traits_ls.txt | \
while read trait; do

###
outdir2=./2_impute.outs/${trait}
zcat ${outdir2}/*_miss.txt.gz |grep -v id_b38 | bgzip > ${outdir2}/all_gwas_impute.txt.gz  

echo ${trait}

done 
