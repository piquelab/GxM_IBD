#!/bin/bash/

cd $PWD

cat traits_ls.txt | \
while read trait; do
###
  dir=./1_missing.outs/${trait}
  if [ ! -d ${dir} ]; then
     mkdir -p ${dir}
  fi 
  # rm ${dir}/zzz_splitSNP*

  cat ${dir}_missing.txt |sed '1d' |split -l 20000 -d -a 4 - ${dir}/zzz_splitSNP
  sleep 10;
  ls ${dir} | grep zzz_splitSNP > ${trait}_missing_snpFile.txt 
  echo ${trait} 
##
done
