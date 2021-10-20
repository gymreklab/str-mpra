#!/bin/bash

INPUT=$1

# Extract the read ID and read sequence, then merge to one line
# Require bases 21-33 match expected filler sequence
# Output the read ID and barcode
zcat $INPUT | awk '(NR%4==1 || NR%4==2)' | \
    awk 'NR % 2 == 1 { o=$0 ; next } { print o "\t" $0 }' | \
    awk -F"\t" '{print $1 "\t" substr($2,1,20) "\t" substr($2,21,33)}' | \
    awk -F"\t" '($3=="TCTAGAGGTTCGTCGACGCGATCGACAGAGACC")' | \
    cut -f 1-2

exit 0
