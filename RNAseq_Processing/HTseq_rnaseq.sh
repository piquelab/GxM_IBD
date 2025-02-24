#!/bin/bash

#SBATCH -q express -p erprp
#SBATCH --mem=100G
#SBATCH --time=1-01:00:00
#SBATCH -N 1
#SBATCH -n 8




module unload python
module load anaconda3.python
source activate htseq


gtffile=/wsu/home/groups/piquelab/pbmc_handls/BAMs/Reference/Homo_sapiens.GRCh38.103.chr.gtf.gz

htseq-count ../bams/${var}_clean.bam  $gtffile --stranded=reverse -f bam  > counts2/${var}.cnts

echo ${var} >> Finished.txt