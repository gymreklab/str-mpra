#!/bin/bash

PREFIX1=new-mid_S7 #$1
PREFIX2=$(echo $PREFIX1 | cut -d '_' -f1)

#BAM=/storage/mlamkin/projects/eSTR-MPRA-analysis/amplification_validation/bams/oligos/default_bwa_params/${PREFIX1}.sorted.bam
BAM=/storage/mlamkin/projects/eSTR-MPRA-analysis/amplification_validation/bams/round1/updated_bwa_params/${PREFIX1}.sorted.bam
#REFFA=/storage/resources/dbase/human/hg19/Illumina/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
REFFA=/storage/resources/dbase/human/hg19/hg19.fa

# Add read group
java -jar /storage/resources/source/picard.jar AddOrReplaceReadGroups \
    I=${BAM} \
    O=${PREFIX2}.rg.bam \
    RGID=${PREFIX2} \
    RGLB=${PREFIX2} \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=${PREFIX2}

samtools index ${PREFIX2}.rg.bam

# Left align
java -jar /storage/resources/source/GenomeAnalysisTK.jar \
    -R ${REFFA} \
    -T LeftAlignIndels \
    -I ${PREFIX2}.rg.bam \
    -o ${PREFIX2}.LA.bam

# Annotate and filter for visualization
./annotate_bam.py \
    ${PREFIX2}.LA.bam ${PREFIX2}.annot.bam
samtools index ${PREFIX2}.annot.bam
