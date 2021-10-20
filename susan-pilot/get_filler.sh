#!/bin/bash

INPUT=$1

zcat $INPUT | awk '(NR%4==1 || NR%4==2)' | \
    awk 'NR % 2 == 1 { o=$0 ; next } { print o "\t" $0 }' | \
    awk -F"\t" '{print $1 "\t" substr($2,21,length($2))}' | \
    awk -F"\t" '(substr($2,1,10)=="TCTAGAGGTT") {split($2,a,"CAGTG"); print $1 "\t" length(a[1]) "\t" a[1]}' 
