## overall comments
if you want this to be easy to use by other people, you should list what packages need to be installed /
other environment configuration details. For example, I can't run STR-BC_Pre-Processing.py --help 
because I don't have Levenshtein, and I don't know what that is.
if you want other people to be able to use this without talking to you, you should describe what
these steps do. Example: pre-processing - I get this does some QC, but can you add a comment
describing what QC its doing?

# MPRA processing 
- location (snolax): `/storage/MPRA/str-mpra/processing/`

## Association

what's the difference between read1 and read2? they seem to have different defaults, so I assume
there's some difference, but I don't know what

- __Association Read Pre-processing via `STR-BC_Pre-Processing.py`__ # minor: I prefer not to have dashes in file names. Can sometimes require the names to be quoted in shell scripts.
    - usage \
      type `./STR-BC_Pre-Processing.py --help` in terminal for detail, below is sample command for modification
    ```shell
    nohup sh -c -u "/storage/MPRA/str-mpra/processing/STR-BC_Pre-Processing.py \
      --read1 /storage/nextseq_demultiplex/NextSeqDate_Project/MPRA_read1.fastq.gz \
      --read2 /storage/nextseq_demultiplex/NextSeqDate_Project/MPRA_read2.fastq.gz \
      --filetype fastq.gz \
      --bwaref /storage/q5gong/MPRA-Susan/array_probes_human_fullprobe_151bp.fa \
      --outdir /storage/MPRA/hSTR1/Date_Initial_Association/outputs/pre-process/ \
      --R1_match 5 \
      --R1_thres 0 \
      --R2_match 16 \
      --R2_thres 0" > /storage/MPRA/hSTR1/Date_Initial_Association/pre_processing.out 2>&1 &
    ```
- __Barcode-STR Association via `STR-BC_Association_v3.py`__
    - usage \
      type `STR-BC_Association_v3.py --help` in terminal for detail, below is sample command for modification
    ```shell
    nohup sh -c -u "/storage/MPRA/str-mpra/processing/STR-BC_Association_v3.py \
      --bam /storage/MPRA/hSTR1/Date_Initial_Association/outputs/pre-process/filtered.bam \
      --outdir storage/MPRA/hSTR1/Date_Initial_Association/outputs/bc-str/ \
      --len 135 \
      --occurrence 5 \
      --minBarcode 3"> /storage/MPRA/hSTR1/Date_Initial_Association/association.out 2>&1 &
    ```
- note 
    - For pre-processing, `--read1`, `--read2`, `--outdir` should always be modifed
    - For association, `--bam`, `--outdir` should always be modified 
        - the .bam file follow by the flag `--bam` is the output from preprocessing 
    - __make sure to create the output directory beforehand__
        - recommended naming for better organization: \
          __snorlax__
          - in `/storage/MPRA/hSTR1/` create a working directory with the format of `Date_Initial_Association` (e.g. `20220712_QG_Association`)
          - in this working directory, create a `outputs` directory with the `pre-process` directory and `bc-str` directory inside

## Expression 
- via `Expression_Read_Processing.py`
    - usage \
      type `Expression_Read_Processing.py --help` in terminal for detail, below is sample command for modification
    ```shell
    nohup sh -c -u "/storage/MPRA/str-mpra/processing/Expression_Read_Processing.py \
      --seqdir /storage/nextseq_demultiplex/NextSeqDate_Project/ \
      --names /storage/MPRA/hSTR1/Date_Initial_Expression/input_file.txt \
      --filetype fastq.gz \
      --numreplicate 3 \
      --association /storage/MPRA/hSTR1/Date_Initial_Association/outputs/bc-str/association.tsv \
      --DESeq2 /storage/MPRA/str-mpra/processing/MPRA_DESeq2.r \
      --refSTR /storage/MPRA/str-mpra/design/array_probes_split.tsv \
      --tss_str /storage/MPRA/str-mpra/design/tss_str_pairs.tsv \
      --ens_gene /storage/MPRA/str-mpra/design/ens_genes.tab \
      --outdir /storage/MPRA/hSTR1/Date_Initial_Expression/outputs/ \
      --read_length 84 \
      --lev_match 5 \
      --lev_thres 5 \
      --min_count 10\
      --min_barcode 3" > /storage/MPRA/hSTR1/Date_Initial_Expression/expression.out 2>&1 &
    ```
- note 
    - For expression read processing, `--seqdir`, `--names`, `--association`, `--outdir` should always be modifed
        - the .tsv file follow by the flag `--association` is the output from STR-BC association 
    - __make sure to create the output directory beforehand__
        - recommended naming for better organization: \
          __snorlax__
          - in `/storage/MPRA/hSTR1/` create a working directory with the format of `Date_Initial_Expression` (e.g. `20220908_SMB_Expression`)
          - in this working directory, create a `outputs` directory
