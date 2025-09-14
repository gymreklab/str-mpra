#!/bin/bash

set -e

# ============================================================================
# SINGLE SAMPLE STR-MPRA EXPRESSION PIPELINE
# ============================================================================
# Users can modify these variables and rerun for different samples/datasets

# MODIFY THESE FOR YOUR SAMPLE:
DATASET=""                              # Your dataset name
SAMPLE=""                               # Your sample name  
DIR_STORAGE="~"                         # Your data directory
THREADS=8                               # Number of CPU threads
ASSO_TABLE=""                           # Path to your association table

# ============================================================================
# HARDCODED VARIABLES (usually don't need to change)
# ============================================================================

# Adapter and trim params
ADAPT_R1="AGATCGGAAGAGC"
ERROR_R1=0.1
MINL_R1=20
MINL_TRIMMED_R1=10
MAXN_R1=3

# Filtering thresholds
BC_READSNUM_THRESHOLD=10
BC_NUM_THRESHOLD=3

# ============================================================================
# PIPELINE PROCESSING
# ============================================================================

echo "Processing sample: ${DATASET}/${SAMPLE}"

# Create directories
mkdir -p ${DIR_STORAGE}/${DATASET}/{clean,expression/${SAMPLE},expression/results}

# Trim adapters
echo "Trimming adapters..."
cutadapt --cores ${THREADS} --discard-untrimmed --max-n ${MAXN_R1} \
    -a ${ADAPT_R1} --minimum-length=${MINL_R1} -e ${ERROR_R1} -q 10 -O ${MINL_TRIMMED_R1} \
    -o ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1.fq.gz \
    ${DIR_STORAGE}/${DATASET}/fq/${SAMPLE}_R1_001.fastq.gz

# Count barcodes
echo "Counting barcodes..."
zcat ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1.fq.gz | \
    awk '(NR%4==2){print $0}' | \
    sort -S10G --parallel=${THREADS} | \
    uniq -c | \
    awk -v OFS="\t" '{print $2, $1}' | \
    sort -S5G --parallel=${THREADS} -t $'\t' -k1 > \
    ${DIR_STORAGE}/${DATASET}/expression/${SAMPLE}/${SAMPLE}_BC_count.tsv

# Assign barcodes to STRs
echo "Assigning barcodes to STRs..."
join -t $'\t' -1 1 -2 1 ${ASSO_TABLE} \
    ${DIR_STORAGE}/${DATASET}/expression/${SAMPLE}/${SAMPLE}_BC_count.tsv > \
    ${DIR_STORAGE}/${DATASET}/expression/${SAMPLE}/${SAMPLE}_BC_str.tsv

# Filter by read number
echo "Filtering by read number..."
awk -v t=${BC_READSNUM_THRESHOLD} '($3>=t){print $0}' \
    ${DIR_STORAGE}/${DATASET}/expression/${SAMPLE}/${SAMPLE}_BC_str.tsv > \
    ${DIR_STORAGE}/${DATASET}/expression/${SAMPLE}/${SAMPLE}_BC_str_filtered.tsv

# Count STRs and barcode numbers
echo "Counting STRs and barcode numbers..."
awk -v FS="\t|_" -v OFS="_" '{print $2, $3, $4, $5, $6 "\t" $NF}' \
    ${DIR_STORAGE}/${DATASET}/expression/${SAMPLE}/${SAMPLE}_BC_str_filtered.tsv | \
    sort -t $'\t' -k1,1 | \
    awk -v OFS="\t" 'BEGIN{BC_count=1} { 
        if (n == $1) { 
            f += $2; BC_count +=1 
        } else { 
            if (n) print n, f, BC_count; 
            n = $1; f = $2; BC_count=1 
        } 
    } END { 
        if (n) print n, f, BC_count 
    }' > ${DIR_STORAGE}/${DATASET}/expression/${SAMPLE}/${SAMPLE}_str_count_BC_num.tsv

# Final filtering by barcode number
echo "Final filtering by barcode number..."
awk -v OFS="\t" -v t=${BC_NUM_THRESHOLD} '($3 >= t) {print $1, $2}' \
    ${DIR_STORAGE}/${DATASET}/expression/${SAMPLE}/${SAMPLE}_str_count_BC_num.tsv > \
    ${DIR_STORAGE}/${DATASET}/expression/results/${SAMPLE}_str_count.tsv

echo "Completed processing sample: ${DATASET}/${SAMPLE}"