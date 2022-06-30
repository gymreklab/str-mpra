#!/usr/bin/env python3
"""
script for processing reads for STR-BC association
"""

# Imports 
import os
import sys
import copy
import gzip
import argparse
import subprocess
import Levenshtein

import numpy as np
import pandas as pd

# Helper functions
def parse_reads_to_df(file_type, read_type, path):
    
    """
    to load reads
    parameter: 
        1. file type of reads - read.fasta.gz, read.fasta, read.fa,
                                read.fastq.gz, read.fastq, read.fq
        2. read type - R1 or R2
        3. path to read file 
    output:
        a data frame that parse the reads 
    """ 
    
    reads = {"read_id":[],
             "sequence":[],
             "quality":[]}

    if file_type in ["fasta.gz", "fastq.gz", "fasta", "fastq", "fa", "fq"]:
        if ".gz" in file_type:
            with gzip.open(path, "rt") as file:
                df = pd.read_csv(file, sep='\n', header=None)
                df = pd.DataFrame(df.values.reshape(-1, 4))
                df.columns=['read_id', read_type, '+', 'qual']
            
            file.close()
            
            return df
        else:
            df = pd.read_csv(path, sep='\n', header=None)
            df = pd.DataFrame(df.values.reshape(-1, 4))
            df.columns=['read_id', read_type, '+', 'qual']
            
            return df
    else:
        raise ValueError("File type is not fasta.gz, fastq.gz, fasta, fa, fastq, or fq")
        

def read_filtering (read_df, read_type, 
                    levenshtein_matching_length=5,
                    levenshtein_threshold=0):

    """
    to filter reads
    parameter:
        1. read_df - data frame that parse the read 
        2. read type - R1 or R2
        3. levenshetein_matching_length - the length of the sequence that will be 
                                          used to calculate levenshetein distancing 
        4. levenshetein_threshold - maximum levenshetein distance the read could have 
    output:
        a data frame with filtered reads
    """    
    
    txt = ("{filt} reads are filtered from {total} reads of {R}, " + 
           "with a levenshtein matching length of {lev_len} and " +
           "threshold of {lev_thres}, " + 
           "resulting in {remain} reads")
    
    reads = copy.deepcopy(read_df)
    reads["length"] = reads[read_type].str.len()
    Read_id = reads["read_id"].str.split(" ", expand=True)
    reads["ID"] = Read_id[0]
    
    lev_scores = []
    # control sequence that ideally should be identical in each read
    # qc_seq
    qc_seq_R1 = "TCTAGAGGTTCGTCGACGCGATCGACAGAGACC" # gibson homology, 33 bp
    qc_seq_R2 = "AACTGGCCGCTTGACG" # adapter sequence, with an extra A from primer, 16 bp
    
    if read_type == "R1":
        qc_seq = qc_seq_R1   
        ref_seq = qc_seq[:levenshtein_matching_length]
        
        matches = []
        
        for seq in reads[read_type]:
            seq_check = seq[20: 20 + levenshtein_matching_length]
            
            # calculate levenshtein distancing
            lev_scores.append(Levenshtein.distance(ref_seq, seq_check))
            
            # check for 2bp exact match 
            if seq_check[:2] == ref_seq[:2]:
                matches.append("matched")
            else:
                matches.append("unmatched")
        
        reads["lev_score"] = lev_scores
        reads["2bp match"] = matches
        total_read = len(reads)
        
        reads = reads[reads["lev_score"] <= levenshtein_threshold]
        reads = reads[reads["2bp match"] == "matched"]
        
        print(txt.format(filt = total_read - len(reads),
                         total = total_read, 
                         R = read_type,
                         lev_len = levenshtein_matching_length,
                         lev_thres = levenshtein_threshold,
                         remain = len(reads)))
        return reads
        
    elif read_type == "R2":
        qc_seq = qc_seq_R2
        ref_seq = qc_seq[:levenshtein_matching_length]
        
        for seq in reads[read_type]:
            seq_check = seq[: levenshtein_matching_length]          
            lev_scores.append(Levenshtein.distance(ref_seq, seq_check))
            
        reads["lev_score"] = lev_scores
        total_read = len(reads)
        
        reads = reads[reads["lev_score"] <= levenshtein_threshold]
        
        print(txt.format(filt = total_read - len(reads),
                         total = total_read, 
                         R = read_type,
                         lev_len = levenshtein_matching_length,
                         lev_thres = levenshtein_threshold,
                         remain = len(reads)))
        return reads
        
    else:
        raise ValueError("Unrecognized read type, please specify R1 or R2")
        
def output_filt_reads (read_df, out_dir, read_type, fname):
    
    """
    output filtered reads
    parameter:
        1. read_df - data frame with filtered reads
        2. out_dir - directory to output the filtered read file 
        3. read_type - R1 or R2
        4. fname - filename 
    output:
        read file with filtered reads 
    """    
    
    if ".gz" in fname:
        file = gzip.open(out_dir + fname, "wb")

        for read_id, seq, sign, qual in zip(read_df["read_id"],
                                            read_df[read_type], 
                                            read_df["+"],
                                            read_df["qual"]):
            file.write((read_id + "\n" + seq + "\n" + 
                        sign + "\n" + qual + "\n").encode())

        file.close()

    else:
        file = open(out_dir + fname, "w")

        for read_id, seq, sign, qual in zip(read_df["read_id"],
                                            read_df[read_type], 
                                            read_df["+"],
                                            read_df["qual"]):
            file.write(read_id + "\n" + seq + "\n" + 
                       sign + "\n" + qual + "\n")

        file.close()

def bwamem_alignment (ref_path, read2_path, out_path):

    """
    perform alignment using bwamem
    parameter:
        1. ref_path - path to reference file 
        2. read2_path - path to read 2 file with STR info
        3. out_path - output bam result
    output:
        bam file of aligned read 2 
    """    
    
    # check if index is already created 
    index_exist = True
    for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
        if not os.path.exists(ref_path + suffix):
            index_exist = False
    
    if not index_exist:
        subprocess.call(("bwa index {bwa_ref}")
                        .format(bwa_ref = ref_path),
                        shell=True)
    
    with open(out_path, "w") as bam_file:
        subprocess.call(("bwa mem -t 6 -L 100 -k 8 -O 5 {bwa_ref} {reads} |\
                         /usr/local/bin/samtools view -bS")
                        .format(bwa_ref = ref_path, 
                                reads = read2_path), 
                        shell=True, stdout=bam_file)  
        
def bam_to_tsv (bam_path, out_path):
        
    """
    processed the bam file
    parameter:
        1. bam_path - path to bam file 
        2. out_path - path to output processed result 
    output:
        processed bam file, stored in tsv format 
    """    
    
    with open(out_path, "w") as tsv_file:
        subprocess.call(("/usr/local/bin/samtools view {bam} | cut -f1,3,6,10")
                        .format(bam = bam_path), 
                        shell=True, stdout=tsv_file)  
        
def getargs():
    parser = argparse.ArgumentParser()
    
    # input file & output directory 
    inout_group = parser.add_argument_group("Input/Output")
    inout_group.add_argument("--read1", help="Path to read1 file", type=str,
                             required=True)
    inout_group.add_argument("--read2", help="Path to read2 file", type=str,
                             required=True)
    inout_group.add_argument("--filetype", help="File type of reads", type=str, 
                             choices=["fasta.gz", "fastq.gz", "fasta", "fastq", "fa", "fq"],
                             required=True)
    inout_group.add_argument("--bwaref", help="Path to ref.fa", type=str,
                             required=True)
    inout_group.add_argument("--outdir", help="Path to output directory", type=str,
                             required=True)
    
    # levenshtein filter 
    filter_group = parser.add_argument_group("Levenstein filtering threshold")
    filter_group.add_argument("--R1_match",
                              help="Length of read 1 that will be use for filter, default is 5",
                              type=int, default=5)
    filter_group.add_argument("--R1_thres",
                              help="Maximum levenshtein score required to kept the read, default is 0",
                              type=int, default=0)
    filter_group.add_argument("--R2_match",
                              help="Length of read 2 that will be use for filter, default is 16",
                              type=int, default=16)
    filter_group.add_argument("--R2_thres",
                              help="Maximum levenshtein score required to kept the read, default is 0",
                              type=int, default=0)
    
    # get argument
    args = parser.parse_args()
    
    return args
    
def main(args):
    # input/output
    R1_path = args.read1
    R2_path = args.read2
    file_type = args.filetype
    bwa_ref = args.bwaref
    out_dir = args.outdir

    # filtering related 
    # levenshtein 
    R1_lev_matching_length = args.R1_match
    R1_lev_threshold = args.R1_thres
    R2_lev_matching_length = args.R2_match
    R2_lev_threshold = args.R2_thres
    
    # check file existence 
    if not os.path.exists(R1_path):
        print("Error: %s does not exist"%R1_path)
        return 1
    
    if not os.path.exists(R2_path):
        print("Error: %s does not exist"%R2_path)
        return 1
    
    if not os.path.exists(bwa_ref):
        print("Error: %s does not exist"%bwa_ref)
        return 1
    
    if not os.path.exists(os.path.dirname(os.path.abspath(out_dir))):
        common.WARNING("Error: The output directory {outdir} does not exist"
                       .format(outdir=out_dir))
        return 1

    # process read 1
    R1 = parse_reads_to_df(file_type, "R1", R1_path)
    R1_filt = read_filtering(R1, "R1", 
                             levenshtein_matching_length=R1_lev_matching_length,
                             levenshtein_threshold=R1_lev_threshold)
    R1_suffix = "_lev" + str(R1_lev_matching_length) + "thres" + str(R1_lev_threshold)
    R1_fname = "R1" + R1_suffix

    print("writing filtered read 1 to " + out_dir + R1_fname  + "." + file_type + "\n")
    output_filt_reads(R1_filt, out_dir, "R1", (R1_fname  + "." + file_type))

    # process read 2
    R2 = parse_reads_to_df("fastq.gz", "R2", R2_path)
    R2_filt = read_filtering(R2, "R2", 
                             levenshtein_matching_length=R2_lev_matching_length,
                             levenshtein_threshold=R2_lev_threshold)
    R2_suffix = "_lev" + str(R2_lev_matching_length) + "thres" + str(R2_lev_threshold)
    R2_fname = "R2" + R2_suffix

    print("writing filtered read 2 to " + out_dir + R2_fname  + "." + file_type + "\n")
    output_filt_reads(R2_filt, out_dir, "R2", (R2_fname  + "." + file_type))

    # alignment on read 2
    out_bam = out_dir + R2_fname + ".bam"
    print("start aligning read 2 to " + bwa_ref)
    bwamem_alignment(bwa_ref, (out_dir + R2_fname + "." + file_type), out_bam)
    print("alignment done" + "\n")

    # process bam file 
    out_tsv = out_dir + R2_fname + ".tsv"
    print("start processing aligning read 2")
    bam_to_tsv(out_bam, out_tsv)
    print("read processing done")    
    
    
    
    return 0

def run():
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__":
    run()