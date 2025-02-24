cat /wsu/home/groups/piquelab/IBD_eQTL/FastQTL/Rectum_protein_coding_residuals/GxE_dichotomous_19PCs/1-per-var-list.txt | tr -d '\015' |  while read var; do
    echo "--- Sumbitting -$var---";
    sbatch -q primary -N1-1 -n 4 --mem=30G -t 10000 --job-name=JOB$var --wrap "module load gnu9 R/test_4.3.2; R --vanilla --args $var < test_GxM_presence_absence_19PCs.R  | gzip > output/output.$var.txt.gz"
    sleep 1
done 