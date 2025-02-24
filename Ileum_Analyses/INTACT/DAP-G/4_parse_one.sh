#!/bin/bash
#SBATCH -q primary
##SBATCH -p erprp
#SBATCH --mem=8G
#SBATCH --time=2-01:00:00
#SBATCH -N 1-1
#SBATCH -n 1


indir="./dap-g_outs"
outdir="./dap-g_parse_outs"
if [ ! -d ${outdir} ]; then
   mkdir -p ${outdir}
fi


cat ./geneList/${geneFile} | \
while read ENSG; do
   ##
   if [ -f ${indir}/${ENSG}.out ]; then
   ##
      echo ${ENSG}
      grep '\[' ${indir}/${ENSG}.out > ${outdir}/${ENSG}.model.out
      grep '((' ${indir}/${ENSG}.out | awk -vENSG=${ENSG} '$0=ENSG "\t" $0' > ${outdir}/${ENSG}.SNP.out
      grep '{'  ${indir}/${ENSG}.out > ${outdir}/${ENSG}.cluster.out
   fi
##
done
