# MPRA processing 

TODO:
* put STR-BC plotting in a separate script

## Required pcakges

See `requirements.txt`.

## Pipeline description

The STR MPRA process consists of two sequencing runs:

1. The STR-BC association library is used to determine which STR sequence each barcode is attached to. Read 1 contains barcode information, and read 2 contains the variable region. Preprocessing of this library is performed with script: `STR-BC_PreProcessing.py`. 

2. The EXPR library is used to count observed barcodes in gDNA and cDNA libraries. Preprocessing of this library is performed with the following scripts:

* `Expression_QC_Processing.py`: TODO
* `Expression_Initial_Analysis.py`: TODO

Additional options to each script are described below.

## STR-BC Association

The script `STR-BC_PreProcessing.py` performs read filtering and mapping, and outputs a list of STR-barcode associations. It is recommended to use default (lenient) filtering options at this stage, and further filtering of barcodes can be done at later stages.

Basic usage (using example files in this repo):

```shell
./STR-BC_PreProcessing.py \
  --read1 test_files/test_reads1.fq.gz \
  --read2 test_files/test_reads2.fq.gz \
  --bwaref test_files/array_probes_human_fullprobe_151bp.fa \
  --outprefix test_output/test
```

Required options:
* `--read1 <STR>` and `--read2 <STR>`: give paths to the reads (fastq or fasta format)
* `--bwaref <STR>`: gives the path to the bwa-indexed reference fasta for the array
* `--outprefix <STR>`: prefix to name output files

Additional options for Levenshtein filtering:
* `--R1_match <int>`: Length of read 1 that will be use for levenshtein filter, default is 5
* `--R2_match <int>`: Length of read 2 that will be use for levenshtein filter, default is 16
* `--R1_thres <int>`: Maximum levenshtein score required to kept the read, default is 0
* `--R2_thres <int>`: Maximum levenshtein score required to kept the read, default is 0

Additional arguments for STR-BC filtering:
* `--len <INT>`: Expected read length for CIGAR filtering. By default, infer from raw read 2
* `--occurrence <INT>`: Minimum required occurence for a unique STR-BC pair
* `--minBarcode  <INT>`: Minimum number of unique barcodes required to be associated per STR
* `--downsample <FLOAT>`: Downsample read pairs by this fraction

This script outputs:
* `$outprefix.summary.csv`: contains summary log info
* `$outprefix.filtered.fq.gz`: filtered read 2 fastq file
* `$outprefix.filtered.bam`: aligned read 2
* `$outprefix.filtered.bam.bai`: index of BAM file for aligned read 2
* `$outprefix.raw_association.tsv`: unfiltered list of STR-BC associations
* `$outprefix.association.tsv`: filtered STR-BC associations

Additional arguments for STR-BC filtering:
* `--occurrence <INT>`: Minimum required occurence for a unique STR-BC pair
* `--minBarcode  <INT>`: Minimum number of unique barcodes required to be associated per STR

The following filters are applied:
* Remove barcodes not associated with exactly 1 STR
* Remove STR-BC pairs supported by less than `occurrence` number of reads
* Remove STRs with fewer than `minBarcode` unique barcodes after filtering.

## Expression 
- read qc and processing via `Expression_QC_Processing.py`
    - usage \
      type `Expression_QC_Processing.py --help` in terminal for detail, below is a sample command for modification
    ```shell
    nohup sh -c -u "/storage/q5gong/str-mpra/processing/Expression_QC_Processing.py \
      --seqdir /storage/nextseq_demultiplex/NextSeqDate_Project/ \
      --names /storage/MPRA/hSTR1/Date_Initial_Expression/input_file.txt \
      --filetype fastq.gz \
      --numreplicate 3 \
      --association /storage/MPRA/hSTR1/Date_Initial_Association/outputs/bc-str/association.tsv \
      --outdir /storage/MPRA/hSTR1/Date_Initial_Expression/outputs/ \
      --export_fastq Yes \
      --read_length 101 \
      --lev_match 5 \
      --lev_thres 5 \
      --min_count 10\
      --min_barcode 3" > /storage/MPRA/hSTR1/Date_Initial_Expression/expression.out 2>&1 &
    ```
- note 
    - For expression read processing, `--seqdir`, `--names`, `--association`, `--outdir` should always be modifed
        - the .tsv file follow by the flag `--association` is the output from STR-BC association 
    - __for the output directory__
        - recommended naming for better organization: \
          __snorlax__
          - in `/storage/MPRA/hSTR1/` create a working directory with the format of `Date_Initial_Expression` (e.g. `20220908_SMB_Expression`)
    - __default for `--export_fastq` is No__
          
- initial analysis via `Expression_Initial_Analysis.py`
    - usage \
      type `Expression_Initial_Analysis.py -- help` in terminal for detail, below are sample commands for modification 
      - for hSTR library 
        ```shell
        nohup sh -c -u "/storage/q5gong/str-mpra/processing/Expression_Initial_Analysis.py \
            --mode hSTR \
            --datadir /storage/MPRA/hSTR1/Date_Initial_Expression/outputs/ \
            --numreplicate 3 \
            --DESeq2 /storage/q5gong/str-mpra/processing/MPRA_DESeq2.r \
            --refSTR /storage/q5gong/str-mpra/design/array_probes_split.tsv \
            --outdir /storage/MPRA/hSTR1/Date_Initial_Expression/outputs/ " > /storage/MPRA/hSTR1/Date_Initial_Expression/initial_analysis.out 2>&1 &
        ```
       - for Uber library
        ```shell
        nohup sh -c -u "/storage/q5gong/str-mpra/processing/Expression_Initial_Analysis.py \
            --mode Uber \
            --datadir /storage/MPRA/hSTR1/Date_Initial_Expression/outputs/ \
            --numreplicate 3 \
            --refFrequency /storage/q5gong/str-mpra/design-uber/round1_probes_sequency_reference.csv \
            --outdir /storage/MPRA/hSTR1/Date_Initial_Expression/outputs/ " > /storage/MPRA/hSTR1/Date_Initial_Expression/initial_analysis.out 2>&1 &
        ```
- note
    - For expression initial analysis, `--mode`, `--datadir`, `--numreplicate` and `--outdir` should always be modified