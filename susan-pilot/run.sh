#!/bin/bash

# Get counts for each UMI. Require the rest of the read to be perfect
for fq in $(ls Susan_MPRA/*.fastq.gz)
do
    echo "Counting UMIs for ${fq}..."
    name=$(basename $fq _R1_001.fastq.gz)
    ./get_counts.sh $fq $name.counts
done

# Combine counts + labels to a single tab file
./combine_counts.py \
    $(ls *.counts) > allcounts.tab

# Filter counts
./filter_counts.py allcounts.tab allcounts_filtered.tab
