#!/bin/bash

INPUT=$1

# Merge all fastq lines to one
# make sure adapter seq matches (bp 1-16)
# extract fq with seq/qual scores for bases 17+

zcat $INPUT | awk 'NR % 2 == 1 { o=$0 ; next } { print o "\t" $0 }' | awk 'NR % 2 == 1 { o=$0 ; next } { print o "\t" $0 }' | \
    awk -F"\t" '{print $0 "\t" substr($2,1, 16)}' | \
    awk -F"\t" '($NF=="AACTGGCCGCTTGACG")' | \
    awk -F "\t" '{print $1 "\n" substr($2,17,length($2)) "\n+\n" substr($4,17,length($4))}'

