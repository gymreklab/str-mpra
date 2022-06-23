#!/bin/bash

#for 64bp homology: TCTAGAGGTTCGTCGACGCGATTATTATCATTACTTGTACAGCTCGTCCATGCCGAGAGTGATC
#for 40bp homology: TCTAGAGGTTCGTCGACGCGATTATTATCATTACTTGTAC
#for 20bp homology: TCTAGAGGTTCGTCGACGCG
#for 10bp homology: TCTAGAGGTT
#for  5bp homology: TCTAG

for file in /storage/q5gong/MPRA-Susan/Next-seq/2022-06-07/*.fastq.gz;
do
    prefix=`basename $file .fastq.gz`
    echo "start processing $file"
    
    zcat $file | awk 'NR % 2 == 1 { o=$0 ; next } { print o "\t" $0 }' | awk 'NR % 2 == 1 { o=$0 ; next } { print o "\t" $0 }' | awk -F"\t" '{print $0 "\t" substr($2,21, 10)}' | awk -F"\t" '($NF=="TCTAGAGGTT")' | awk -F "\t" '{print $1 "\t" $2}' > /storage/q5gong/MPRA-Susan/Next-seq/2022-06-07/Processed/${prefix}-proc-10.tsv
    
done