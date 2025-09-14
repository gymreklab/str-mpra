#!/bin/bash

set -e

# ============================================================================
# DUAL ALIGNMENT STR-MPRA PIPELINE
# ============================================================================
# Users can modify these variables and rerun for different samples/datasets

# MODIFY THESE FOR YOUR SAMPLE:
DATASET=""              # Your dataset name
SAMPLE=""               # Your sample name  
MODE="match"            # Alignment mode: "match" or "indel"
IDX_BWA_HG38_FWD=""     # Path to your forward BWA index
IDX_BWA_HG38_REV=""     # Path to your reverse BWA index
DIR_STORAGE=""          # Your data directory
THREADS=8               # Number of CPU threads
REF_FILE=""             # Path to STR_rpt_number_human.csv
TMPDIR="/tmp"           # Path to your tmp folder

# ============================================================================
# HARDCODED VARIABLES (usually don't need to change)
# ============================================================================

READ_THRESHOLD=3
READ_LENGTH=135

# Adapter and trim params
ADAPT_R2="AACTGGCCGCTTGACG"
ADAPT_R1="TCTAG"
ERROR_R2=0.1
ERROR_R1=0
MINL_R2=70
MINL_R1=20
MINL_TRIMMED_R2=9
MAXN_R2=0.1

# Script paths - UPDATE THESE TO YOUR ACTUAL PATHS
CIGAR_FILTER_SCRIPT="CIGAR_filter.py"
DUAL_ALIGNMENT_SCRIPT="dual_alignment_selector_dpSTR.sh"
STUTTER_ERROR_SCRIPT="stutter_error.py"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

check_file() {
    if [[ ! -f "$1" ]]; then
        echo "Error: Required file $1 not found!"
        exit 1
    fi
}

check_script() {
    if [[ ! -f "$1" ]]; then
        echo "Warning: Script $1 not found! Please update the path."
        # Don't exit here as user might have different paths
    fi
}

# ============================================================================
# VALIDATION
# ============================================================================

echo "Validating inputs..."

if [[ -z "$DATASET" || -z "$SAMPLE" ]]; then
    echo "Error: DATASET and SAMPLE must be set!"
    exit 1
fi

check_file "$IDX_BWA_HG38_FWD"
check_file "$IDX_BWA_HG38_REV"
check_file "$REF_FILE"
check_script "$CIGAR_FILTER_SCRIPT"
check_script "$EXTRACT_BAM_SCRIPT"
check_script "$DUAL_ALIGNMENT_SCRIPT"
check_script "$STUTTER_ERROR_SCRIPT"
check_script "$THRESHOLD_KEEP_SCRIPT"

# ============================================================================
# PIPELINE PROCESSING
# ============================================================================

echo "Processing sample: ${DATASET}/${SAMPLE}"

# Create directories
mkdir -p ${DIR_STORAGE}/${DATASET}/{clean,align,strBC/{mergeIDX,stutterError}}

# ============================================================================
# STEP 1: TRIM READS
# ============================================================================

echo "Step 1: Trimming reads..."

# First cutadapt run
cutadapt --cores $THREADS --discard-untrimmed --pair-filter=both --max-n $MAXN_R2 \
    -G X${ADAPT_R2} --minimum-length=$MINL_R2 -e $ERROR_R2 -q 10 -O $MINL_TRIMMED_R2 \
    -o ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1_s1.fq.gz \
    -p ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R2_s1.fq.gz \
    ${DIR_STORAGE}/${DATASET}/fq/${SAMPLE}_R1_001.fastq.gz \
    ${DIR_STORAGE}/${DATASET}/fq/${SAMPLE}_R2_001.fastq.gz

# Second cutadapt run
cutadapt --cores $THREADS --discard-untrimmed --pair-filter=both \
    -a ${ADAPT_R1}$ --minimum-length=$MINL_R1 -e $ERROR_R1 \
    -o ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1.fq.gz \
    -p ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R2.fq.gz \
    ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R1_s1.fq.gz \
    ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}_R2_s1.fq.gz

# ============================================================================
# STEP 2: MERGE BARCODES
# ============================================================================

echo "Step 2: Merging barcodes..."

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

# ============================================================================
# STEP 3: DUAL MAPPING
# ============================================================================

echo "Step 3: Dual mapping (forward and reverse)..."

# Forward mapping
echo "  - Forward mapping..."
bwa mem -M -k 70 -t 8 -O 28 -L 10 $IDX_BWA_HG38_FWD \
    ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}.fq.gz | \
samtools view -@ 8 -bS > ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_fwd.bam

# Reverse mapping
echo "  - Reverse mapping..."
bwa mem -M -k 70 -t 8 -O 28 -L 10 $IDX_BWA_HG38_REV \
    ${DIR_STORAGE}/${DATASET}/clean/${SAMPLE}.fq.gz | \
samtools view -@ 8 -bS > ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_rev.bam

# ============================================================================
# STEP 4: CIGAR FILTERING
# ============================================================================

echo "Step 4: CIGAR filtering..."

# Get sample length from config or use default
SAMPLE_LENGTH=$READ_LENGTH

# Filter forward alignment
python $CIGAR_FILTER_SCRIPT \
    -i ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_fwd.bam \
    -o ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_fwd_filtered.bam \
    -l $SAMPLE_LENGTH -t 21

# Filter reverse alignment
python $CIGAR_FILTER_SCRIPT \
    -i ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_rev.bam \
    -o ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_rev_filtered.bam \
    -l $SAMPLE_LENGTH -t 21

# ============================================================================
# STEP 5: INITIAL ALIGNMENT TABLE
# ============================================================================

echo "Step 5: Creating initial alignment tables..."

# Extract BAM info for forward
bash $EXTRACT_BAM_SCRIPT \
    --length $READ_LENGTH \
    --input ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_fwd_filtered.bam \
    --output ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_fwd.csv

# Extract BAM info for reverse
bash $EXTRACT_BAM_SCRIPT \
    --length $READ_LENGTH \
    --input ${DIR_STORAGE}/${DATASET}/align/${SAMPLE}_rev_filtered.bam \
    --output ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_rev.csv

# ============================================================================
# STEP 6: SORT ALIGNMENT TABLES
# ============================================================================

echo "Step 6: Sorting alignment tables..."

sort -k1,1 ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_fwd.csv > \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_fwd_sorted.csv

sort -k1,1 ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_rev.csv > \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_rev_sorted.csv

# ============================================================================
# STEP 7: MERGE ALIGNMENT TABLES
# ============================================================================

echo "Step 7: Merging alignment tables..."

join -j 1 -a 1 -a 2 -t"," -e"NA" -o 1.1,1.2,1.3,1.4,1.5,1.6,2.1,2.2,2.3,2.4,2.5,2.6 \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_fwd_sorted.csv \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_rev_sorted.csv > \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_merged.csv

# ============================================================================
# STEP 8: DUAL ALIGNMENT FILTER
# ============================================================================

echo "Step 8: Dual alignment filtering..."

bash $DUAL_ALIGNMENT_SCRIPT \
    -m $MODE -l $READ_LENGTH \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_merged.csv \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_${MODE}_fwd_rev_agree.tsv \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_${MODE}_fwd_rev_disagree.tsv

# Process both match and indel modes if needed
for mode in "match" "indel"; do
    if [[ -f $DUAL_ALIGNMENT_SCRIPT ]]; then
        bash $DUAL_ALIGNMENT_SCRIPT \
            -m $mode -l $READ_LENGTH \
            ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_merged.csv \
            ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_${mode}_fwd_rev_agree.tsv \
            ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_${mode}_fwd_rev_disagree.tsv
    fi
done

# ============================================================================
# STEP 9: RAW ASSOCIATION TABLE
# ============================================================================

echo "Step 9: Creating raw association table..."

# Combine all valid alignments
cat ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_match_fwd_rev_agree.tsv \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_match_fwd_rev_disagree.tsv \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_indel_fwd_rev_agree.tsv \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_indel_fwd_rev_disagree.tsv > \
    ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_all_valid_alignment.tsv

# Create association table
cut -f1,3 ${DIR_STORAGE}/${DATASET}/strBC/mergeIDX/${SAMPLE}_all_valid_alignment.tsv | \
awk -v FS=":|\t|," '
    NR==FNR { nrpt[$1] = $2; next } 
    {
        if ($(NF-1) ~ /^chr/) 
        { print $(NF-1) ":" $NF, $(NF-2) } 
        else {
            split($(NF), s, "_");
            print "chrN:start-end_Human-STR-" s[3] "_" s[4] "-original_" nrpt[$NF] "_" s[4] "_0", $(NF-1);
        }
    }
' $REF_FILE - > ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_asso.tsv

# Create unique association table
sort ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_asso.tsv | \
uniq -c | sort -nr | awk -v OFS="\t" '{print $1,$2,$3}' > \
${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_uniq_asso.tsv

# ============================================================================
# STEP 10: STUTTER ERROR CORRECTION BY REPEAT LENGTH
# ============================================================================

echo "Step 10: Stutter error correction by repeat length..."

for i in {1..6}; do
    echo "  - Processing length $i..."
    
    # Filter by length
    awk -F "_" -v l=$i '{ 
        if ($5 == "ref") { 
            if (length($3) == l && $4!=0) { print $0 } 
        } else if ($5 !~ /random/) { 
            if (length($5)==l && $4!=0) { print $0 } 
        } 
    }' ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_uniq_asso.tsv > \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_${i}
    
    # Apply stutter correction
    python $STUTTER_ERROR_SCRIPT \
        -i ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_${i} \
        -o ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_${i}_stutter_corrected \
        -u 0.05 -d 0.05 -m uber \
        --dup_bcs_output ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_${i}_invalid_BC.tsv
done

# ============================================================================
# STEP 11: SPLIT RANDOM AND ZERO MOTIF
# ============================================================================

echo "Step 11: Processing random and zero motif sequences..."

awk -F "[\t_]" '(($0 ~ /random/) || $5==0 ){print $0}' \
    ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_filtered_uniq_asso.tsv | 
awk -v OFS="\t" '{print $1, $3, $2}' | sed 's/random-matchedGC/randomMatchedGC/g' | 
awk -F "[-\t:_]" -v OFS="\t" '
{
    if ($0 ~ /original/){
        print $1, $2, $6 "_" $7 "_" $8 "-" $9 "-" $10 "-" $13 "_" $12 "_p" $11
    }
    else {
        print $1, $2, $6 "_" $7 "_" $8 "-" $9 "-" $11 "-" $12 "_" $11 "_p" $10
    }
}' > ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_random_or_zero

# ============================================================================
# STEP 12: MERGE STUTTER CORRECTED TABLE
# ============================================================================

echo "Step 12: Merging stutter corrected tables..."

# Merge all stutter corrected files with random motifs
cat ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_*_stutter_corrected \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_random_or_zero > \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected.tsv

# Remove all duplicates
echo "Remove all dups"
cut -f2 ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected.tsv | \
sort | uniq -c | sort -k1nr | awk '($1>1){print $2}' > \
${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_dup_bc.txt

awk -v FS="\t" -v OFS="\t" '(NR==FNR){h[$1] = 1; next} !($2 in h){print $1,$2,$3}' \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_dup_bc.txt \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected.tsv > \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_noDupBC.tsv

# ============================================================================
# STEP 13: REFORMAT CORRECTED MERGED FILE
# ============================================================================

echo "Step 13: Reformatting corrected merged file..."

awk -v FS='[\t\r]' -v OFS='\t' '{print $2, $3"_"$1}' \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected.tsv > \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_reformatted.tsv

# ============================================================================
# STEP 14: MERGE ALL BCs
# ============================================================================

echo "Step 14: Merging all BCs..."

cat ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_STR_length_*_invalid_BC.tsv \
    ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_reformatted.tsv | \
sort -t$'\t' -k1,1 | tr -d '\r' > \
${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_reformatted_with_dupBC_sorted.tsv

# ============================================================================
# STEP 15: AGGREGATE STRs BY BC
# ============================================================================

echo "Step 15: Aggregating STRs by barcode..."

awk -F'\t' '
{
    if (prev_bc && prev_bc != $1) { # If the barcode changes (and it is not the first line), print the result for the previous barcode
        print prev_bc "\t" a[prev_bc]  # Print the previous barcode and its concatenated second-column values
    }

    a[$1] = (a[$1] ? a[$1] "," $2 : $2) # Concatenate second column values for the current barcode. If the barcode has appeared before, append the second column value, otherwise initialize the list
    prev_bc = $1
}
END {
    if (prev_bc) {
        print prev_bc "\t" a[prev_bc]
    }
}' ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_reformatted_with_dupBC_sorted.tsv > \
${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_results.tsv

# Separate unique and duplicate BCs
grep , ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_results.tsv > \
${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_dupBC.tsv

grep -v , ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_results.tsv > \
${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_uniqBC.tsv

# ============================================================================
# STEP 16: DUPLICATE BC FILTER
# ============================================================================

echo "Step 16: Filtering duplicate barcodes..."

bash $THRESHOLD_KEEP_SCRIPT \
    -i ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_dupBC.tsv \
    -s ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_dupBC_keep.tsv \
    -f ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_dupBC_remove.tsv \
    -p 99

# ============================================================================
# STEP 17: FINAL ASSOCIATION TABLE WITH DUPLICATE BC
# ============================================================================

echo "Step 17: Creating final association table with duplicate BC handling..."

cut -f1,2 ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_dupBC_keep.tsv | \
cat - ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_stutter_corrected_barcode_aggregated_uniqBC.tsv | \
awk -F'\t' '{ 
    split($2, arr, "_"); 
    value = arr[length(arr)]; 
    sub("_"value"$", "", $2); 
    print value "\t" $1 "\t" $2 
}' > ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_keep_dominateBC_asso_noDupBC.tsv

# ============================================================================
# STEP 18: FINAL BC TABLE NO DUPLICATE BCs
# ============================================================================

echo "Step 18: Creating final BC table without duplicate BCs..."

sed 's/\r//g' ${DIR_STORAGE}/${DATASET}/strBC/stutterError/${SAMPLE}_asso_noDupBC.tsv | \
awk -v OFS="\t" -v t=$READ_THRESHOLD '($1>=t) {print $2, $3 "_" $1}' | \
sort -S5G -t $'\t' -k1 > ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_BC_STRs.tsv

# ============================================================================
# CLEANUP (optional)
# ============================================================================

echo "Pipeline completed successfully!"
echo "Final output: ${DIR_STORAGE}/${DATASET}/strBC/${SAMPLE}_BC_STRs.tsv"
