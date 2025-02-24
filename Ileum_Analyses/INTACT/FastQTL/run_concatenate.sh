#!/bin/bash



cd fastQTL.output
tissue=IBD_ileum
zcat ${tissue}_nominals_chunk*.txt.gz  | gzip -c > ${tissue}_nominals_all_chunks.txt.gz  
