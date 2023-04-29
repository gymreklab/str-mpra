# MPRA processing 
- location (snolax): `/storage/MPRA/str-mpra/processing/`

## Package needed 
- `os`
- `sys`
- `copy`
- `gzip`
- `numpy`
- `Cigar`
- `pandas`
- `seaborn`
- `argparse`
- `subprocess`
- `matplotlib`
- `Levenshtein`
- `scipy.stats`

## Association
- __Association Read Pre-processing via `STR-BC_Pre-Processing.py`__
    - usage \
      type `STR-BC_Pre-Processing.py --help` in terminal for detail, below is a sample command for modification
    ```shell
    nohup sh -c -u "/storage/q5gong/str-mpra/processing/STR-BC_Pre-Processing.py \
      --read1 /storage/nextseq_demultiplex/NextSeqDate_Project/MPRA_read1.fastq.gz \
      --read2 /storage/nextseq_demultiplex/NextSeqDate_Project/MPRA_read2.fastq.gz \
      --filetype fastq.gz \
      --bwaref /storage/q5gong/MPRA-Susan/array_probes_human_fullprobe_151bp.fa \
      --outdir /storage/MPRA/hSTR1/Date_Initial_Association/outputs/ \
      --R1_match 5 \
      --R1_thres 0 \
      --R2_match 16 \
      --R2_thres 0" > /storage/MPRA/hSTR1/Date_Initial_Association/pre_processing.out 2>&1 &
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
            --outdir /storage/MPRA/hSTR1/Date_Initial_Expression/outputs/ " > /storage/MPRA/hSTR1/Date_Initial_Expression/initial_analysis.out 2>&1 &
        ```
- note
    - For expression initial analysis, `--mode`, `--datadir`, `--numreplicate` and `--outdir` should always be modified