# MPRA processing 
- location (snolax): `/storage/MPRA/str-mpra/processing/`

## Required pcakges

See `requirements.txt`.

## Pipeline description

The STR MPRA process consists of two sequencing runs:

1. The STR-BC association library is used to determine which STR sequence each barcode is attached to. Read 1 contains barcode information, and read 2 contains the variable region. Preprocessing of this library is performed with the following scripts:

* `STR-BC_Pre-Processing.py`: performs mapping and filtering of raw reads
* `STR-BC_Association_v3.py`: obtains the STR-BC association file `associations.tsv`

2. The EXPR library is used to count observed barcodes in gDNA and cDNA libraries. Preprocessing of this library is performed with the following scripts:

* `Expression_QC_Processing.py`: TODO
* `Expression_Initial_Analysis.py`: TODO

Additional options to each script are described below.

## STR-BC Association

### Step 1: Mapping and filtering raw reads

Basic usage (using example files in this repo):

```shell
./STR-BC_Pre-Processing.py \
  --read1 test_files/test_reads1.fq.gz \
  --read2 test_files/test_reads2.fq.gz \
  --bwaref test_files/array_probes_human_fullprobe_151bp.fa \
  --outdir test_output/
```

Required options:
* `--read1` and `--read2`: give paths to the reads (fastq or fasta format)
* `--bwaref`: gives the path to the bwa-indexed reference fasta for the array
* `--outdir`: gives the output directory to store results in

Additional options for Levenshtein filtering:

* `--R1_match <int>`: Length of read 1 that will be use for levenshtein filter, default is 5
* `--R2_match <int>`: Length of read 2 that will be use for levenshtein filter, default is 16
* `--R1_thres <int>`: Maximum levenshtein score required to kept the read, default is 0
* `--R2_thres <int>`: Maximum levenshtein score required to kept the read, default is 0

This script outputs:
* `summary.csv` contains the total number of read pairs considered, the number of read pairs remaining after filtering, and the perfect of read pairs remaining after filtering
* `filtered.fq.gz`: filtered read 2 fastq file
* `filtered.bam`: aligned read 2

### Step 2: Obtain STR_BC associations

Basic usage (using example files in this repo):

```shell


```

- __Barcode-STR Association via `STR-BC_Association_v3.py`__
    - usage \
      type `STR-BC_Association_v3.py --help` in terminal for detail, below is a sample command for modification
    ```shell
    nohup sh -c -u "/storage/q5gong/str-mpra/processing/STR-BC_Association_v3.py \
      --bam /storage/MPRA/hSTR1/Date_Initial_Association/outputs/pre-process/filtered.bam \
      --outdir storage/MPRA/hSTR1/Date_Initial_Association/outputs/ \
      --len 135 \
      --occurrence 5 \
      --minBarcode 3"> /storage/MPRA/hSTR1/Date_Initial_Association/association.out 2>&1 &
    ```
- note 
    - For pre-processing, `--read1`, `--read2`, `--outdir` should always be modifed
    - For association, `--bam`, `--outdir` should always be modified 
        - the .bam file follow by the flag `--bam` is the output from preprocessing 
    - __for the output directory__
        - recommended naming for better organization: \
          __snorlax__
          - in `/storage/MPRA/hSTR1/` create a working directory with the format of `Date_Initial_Association` (e.g. `20220712_QG_Association`)
          - in this working directory, create a `outputs` directory with the `pre-process` directory and `bc-str` directory inside

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