#!/bin/bash

INFILE=$1
OUTFILE=$2

zcat $INFILE | awk '(NR%4==2)' | grep "TCTAGAGGTTCGTCGACGCGATTATTATCATTACTTGTACAGCTCGTCCATGCCGAGAGTGATC" | \
    cut -c 1-20 | sort | uniq -c | awk '{print $2 "\t" $1}' > $OUTFILE
