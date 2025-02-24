#!/bin/bash/

cd $PWD

cat traits_ls.txt | \
while read trait; do
   ##
   cat ${trait}_missing_snpFile.txt | \
   while read snpfn; do

   sbatch -q express -p erprp --mem=20G --time=3-00:00:00 -N 1-1 -n 1 --job-name=impute_${trait}_${snpfn} --output=slurm_impute_${trait}_${snpfn}.output --wrap "
   module load R; 
   R CMD BATCH --no-save --no-restore '--args ${trait} ${snpfn}' 2_impute.R 2_impute_${trait}_${snpfn}.Rout"
   echo ${trait} ${snpfn}
   sleep 1;
   done 
   sleep 2;
done 
