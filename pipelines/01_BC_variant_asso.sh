#!/bin/bash

set -e

# ============================================================================
# SINGLE SAMPLE hSTR PIPELINE
# ============================================================================
# Users can modify these variables and rerun for different samples/datasets

# MODIFY THESE FOR YOUR SAMPLE:
DATASET=""              # Your dataset name
SAMPLE=""                  # Your sample name  
IDX_BWA_HG38="strs_probe_hg38_100k.fa"  # Path to your BWA index
DIR_STORAGE="~"                          # Your data directory
THREADS=8                                # Number of CPU threads
CIGAR_FILTER="CIGAR_filter.py"          # Path to CIGAR filter script
STUTTER_ERROR="stutter_error.py"    # Path to stutter error script
TMPDIR="." # Path to your tmp folder

# ============================================================================
# HARDCODED VARIABLES (usually don't need to change)
# ============================================================================

READ_THRESHOLD=3
SAMPLE_LENGTH=135

# Adapter and trim params
ADAPT_R2="AACTGGCCGCTTGACG"
ADAPT_R1="TCTAG"
ERROR_R2=0.1
ERROR_R1=0
MINL_R2=70
MINL_R1=20
MINL_TRIMMED_R2=9
MAXN_R2=0.1

# ============================================================================
# PIPELINE PROCESSING
# ============================================================================

echo "Processing sample: ${DATASET}/${SAMPLE}"

# Create directories
mkdir -p ${DIR_STORAGE}/${DATASET}/{clean,align,strBC/stutterError}

# Trim reads
echo "Trimming reads..."
cutadapt --cores $THREADS --discard-untrimmed --pair-filter=both --max-n $MAXN_R2 \
    -G X${ADAPT_R2} --minimum-length=$MINL_R2 -e $ERROR_R2 -q 10 -O $MINL_TRIMMED_R2 \
    -o ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1_s1.fq.gz \
    -p ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R2_s1.fq.gz \
    ${DIR_STORAGE}/${DATASET}/fq/${SAMPLE}_R1_001.fastq.gz \
    ${DIR_STORAGE}/${DATASET}/fq/${SAMPLE}_R2_001.fastq.gz

cutadapt --cores $THREADS --discard-untrimmed --pair-filter=both \
    -a ${ADAPT_R1}$ --minimum-length=$MINL_R1 -e $ERROR_R1 \
    -o ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1.fq.gz \
    -p ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R2.fq.gz \
    ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1_s1.fq.gz \
    ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R2_s1.fq.gz

# Merge barcodes
echo "Merging barcodes..."
awk '(NR==FNR) {
    if (NR % 4 == 1) {
        read = $1
    }
    if (NR % 4 == 2) {
        BC[read] = $0
    }
}
(NR!=FNR) {
    if (NR % 4 == 1) {
        print $1 ":BC:" BC[$1] " " $2
    } else {
        print $0
    }
}' <(gzip -dc ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1.fq.gz) \
  <(gzip -dc ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R2.fq.gz) | \
gzip > ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}.fq.gz

# Mapping
echo "Mapping reads..."
bwa mem -M -t $THREADS $IDX_BWA_HG38 ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}.fq.gz | \
samtools view -bS | \
samtools sort -T $TMPDIR -m 4G -@ $THREADS > ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}.bam

# CIGAR filter
echo "CIGAR filtering..."
python $CIGAR_FILTER \
    -i ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}.bam \
    -o ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_filtered.bam \
    -l $SAMPLE_LENGTH -t $THREADS

# Raw association table
echo "Creating association table..."
samtools view ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_filtered.bam | \
cut -f1,3 | awk -v FS=":|\t" '{print $NF, $(NF-1)}' > \
${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_asso.tsv

sort ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_asso.tsv | \
uniq -c | sort -nr | awk -v OFS="\t" '{print $1,$2,$3}' > \
${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_uniq_asso.tsv

rm -f ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_asso.tsv

# Stutter error correction by length
echo "Stutter error correction..."
for i in {1..6}; do
    awk -F "_" -v l=$i 'length($4) == l' \
        ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_uniq_asso.tsv > \
        ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_${i}
    
    python $STUTTER_ERROR \
        -i ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_${i} \
        -o ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_${i}_stutter_corrected \
        -u 0.05 -d 0.05
done

# Merge stutter corrected
echo "Merging stutter corrected data..."
cat ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_*_stutter_corrected > \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected.tsv

cut -f2 ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected.tsv | \
sort | uniq -c | sort -k1nr | awk '($1>1){print $2}' > \
${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_dup_bc.txt

awk -v FS="\t" -v OFS="\t" '(NR==FNR){h[$1] = 1; next} !($2 in h){print $1,$2,$3}' \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_dup_bc.txt \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected.tsv > \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_noDupBC.tsv

# Final BC table
echo "Creating final BC table..."
sed 's/\r//g' ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_noDupBC.tsv | \
awk -v OFS="\t" -v t=$READ_THRESHOLD '($1>=t) {match($3, /^([^_]+_[^_]+_[^_]+_[^_]+_[^_]+)/, m); print $2, m[1]}' |\
sort -S5G -t $'\t' -k1 > ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_BC_STRs.tsv

echo "Completed processing sample: ${DATASET}/${SAMPLE}"