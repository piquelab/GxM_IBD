#!/bin/bash



cd fastQTL.output

zcat IBD_Rectum_nominals_chunk*.txt.gz | gzip -c > IBD_Rectum_nominals_all_chunks.txt.gz  
