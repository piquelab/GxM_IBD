Folder contains bam files that have been aligned with HiSAT2. Quality control and de-duplication accomplished using samtools.


Commands:

## ls /wsu/home/fx/fx78/fx7820/rprscratch/Francesca/PBMC/PBP3/fastqs/PBP3-HT* |cut -d/ -f12 | cut -d_ -f 1|sort|uniq > names.txt #Creates file of names of samples to analyze. Used in qsub command.
#for i in `cat names.txt`; do while :; do if [ `qme | wc -l` -lt 120 ]; then qsub -v Sample=${i} align-hisat.sh; break; else sleep 1; fi; done; done

# Replace _ with -

#Making counts.csv
for i in `ls *sorted_count.txt`; do echo "$i"| cut -d_ -f1 >m ; cat m| tr "\n" " "; cat $i; done > sorted.txt
for i in `ls *quality_count.txt`; do echo "$i"| cut -d_ -f1 >m ; cat m| tr "\n" " "; cat $i; done > quality.txt
for i in `ls *clean_count.txt`; do echo "$i"| cut -d_ -f1 >m ; cat m| tr "\n" " "; cat $i; done > clean.txt

join sorted.txt quality.txt > tmp.txt
join tmp.txt clean.txt > all_counts.txt

rm m
rm tmp.txt
rm sorted.txt
rm quality.txt
rm clean.txt

#To make list of total number of reads processed by HiSAT
#ls BCol1* | cut -d_ -f1 | sort | uniq > names.txt
for i in `cat names2.txt`; do echo -n "$i " >> total_reads.txt; head -1 ${i}_aligned.bam.e|cut -d ' ' -f1 >> total_reads.txt;done

mkdir out #for sh.e* and sh.o* files
mv *.out out
mv *.err out


Files:

align.sh: Does alignment, sorting of bam files, quality filtering, de-duplication
names.txt: List of sample names to be used with qsub command
*_aligned.bam: alignment of reads by HiSAT2
*_aligned.bam.e: output of HiSAT2
*_clean.bam: De-duplicated bam file
*_clean.bam.bai: samtools index for clean.bam
*_clean_count.txt: number of reads mapping after removing pcr duplicates
*_sorted_count.txt: number of reads mapped
*_quality.bam: number of reads after quality filtering - removing reads that do not align well to genome
*_quality.bam.bai: index for quality.bam
*_quality_count.txt: number of reads after quality filter (before de-duplication)
*_sorted.bam: sorted alignment bam
*_sorted.bam.bai: index for sorted alignment bam
all_counts.csv: contains all of data from *_count.txt; header=T
total_reads.txt: number of reads processed by HiSAT (before alignment)


Folders:
e_o: contains *.sh.e* and *.sh.o* files
