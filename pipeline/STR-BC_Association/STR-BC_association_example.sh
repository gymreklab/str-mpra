#!/bin/bash

# Read processing
/storage/q5gong/MPRA-Susan/STR-BC_Association_Scripts/STR-BC_read_processing_v2.py \
    --read1 /storage/q5gong/MPRA-Susan/lz0504/lz_S16_L001_R1_001.fastq.gz \
    --read2 /storage/q5gong/MPRA-Susan/lz0504/lz_S16_L001_R2_001.fastq.gz \
    --filetype fastq.gz \
    --bwaref /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/array_probes_human_fullprobe_151bp.fa \
    --outdir /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/read_processing_v2/ \
    --R1_match 5 \
    --R1_thres 0 \
    --R2_match 16 \
    --R2_thres 0 \

# Association
/storage/q5gong/MPRA-Susan/STR-BC_Association_Scripts/STR-BC_Association.py \
    --filtR1 /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/R1_lev5thres0.fastq.gz \
    --tsvR2 /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/R2_lev16thres0.tsv \
    --outdir /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/ \
    --lenR1 151 \
    --lenR2 151 \
    --occurrence 1 \

# example of checking
# processed read 1
echo "processed read1 check"
diff -q <(zcat /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/ipynb_test/R1_lev5thres0.fastq.gz) \
          <(zcat /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/read_processing_v2/R1_lev5thres0.fastq.gz) \
    && echo same || echo not_same
    
# processed read 2
echo "processed read2 check"
diff -q <(zcat /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/ipynb_test/R2_lev16thres0.fastq.gz) \
          <(zcat /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/read_processing_v2/R2_lev16thres0.fastq.gz) \
    && echo same || echo not_same
    
# processed tsv
echo "processed bam check"
diff -q <(cat /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/ipynb_test/R2_lev16thres0.tsv) \
          <(cat /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/read_processing_v2/R2_lev16thres0.tsv) \
    && echo same || echo not_same

# association tsv
echo "association.tsv check"
diff -q <(cat /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/ipynb_test/association.tsv) \
          <(cat /storage/q5gong/MPRA-Susan/pipeline_testing/STR-BC-association/script_test/read_processing_v2/association.tsv) \
    && echo same || echo not_same