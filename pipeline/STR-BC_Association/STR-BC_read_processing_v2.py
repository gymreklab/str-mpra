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
def process(lines=None):
    ks = ['read_id', 'seq', 'optional', 'qual']
    return {k: v for k, v in zip(ks, lines)}

def check_R1 (sequence, levenshtein_matching_length,
              levenshtein_threshold):
    """
    to check read1
    parameter:
        1. sequence
            - the read sequence
        2. levenshetein_matching_length
            - the length of the sequence that will be 
              used to calculate levenshetein distancing 
        3. levenshetein_threshold
            - maximum levenshetein distance the read 
              could have 
    output:
        return true if the sequence pass the check, 
        return false other wise
    """
    
    # the full identical sequence that could be 
    # used for qc
    qc_seq = "TCTAGAGGTTCGTCGACGCGATCGACAGAGACC"
    
    # the reference sequence
    ref_seq = qc_seq[:levenshtein_matching_length]
    
    # the sequence being checked 
    seq_check = sequence[20: 20 + levenshtein_matching_length]
    
    # check for 2bp exact match after barcode
    if seq_check[:2] == ref_seq[:2]:
        
        # calculate levenshtein distancing
        lev_score = Levenshtein.distance(ref_seq, seq_check)
        # check for leveshtein distancing
        if lev_score > levenshtein_threshold:
            # have more mismtach than allowed 
            # fail the check 
            return False
        
    else:
    # does not have the 2bp exact match after barcode
    # fail the check
        return False
    
    # pass the check
    return True

def check_R2 (sequence, levenshtein_matching_length,
              levenshtein_threshold):
    """
    to check read2
    parameter:
        1. sequence
            - the read sequence
        2. levenshetein_matching_length
            - the length of the sequence that will be 
              used to calculate levenshetein distancing 
        3. levenshetein_threshold
            - maximum levenshetein distance the read 
              could have 
    output:
        return true if the sequence pass the check, 
        return false other wise
    """
    
    # the full identical sequence that could be 
    # used for qc
    qc_seq = "AACTGGCCGCTTGACG"
    
    # the reference sequence
    ref_seq = qc_seq[:levenshtein_matching_length]
    
    # the sequence being checked 
    seq_check = sequence[: levenshtein_matching_length]
    
    # calculate levenshtein distancing
    lev_score = Levenshtein.distance(ref_seq, seq_check)
    # check for leveshtein distancing
    if lev_score > levenshtein_threshold:
        # have more mismtach than allowed 
        # fail the check 
        return False
        
    # pass the check
    return True

def filter_read (path, file_type, read_type,
                 out_name, out_dir,
                 lev_matching_length,
                 lev_threshold):
    """
    to filter reads and write the filter read 
    to a new read file of the same file type
    parameter:
        1. R1_path 
            - path to read file
        2. file_type
            - type of read file, fastq.gz or fastq or fq
        3. read_type
            - type of read, R1 or R2
        4. out_name
            - output file name
        5. out_dir
            - output directory 
        6. levenshetein_matching_length
            - the length of the sequence that will be 
              used to calculate levenshetein distancing 
        7. levenshetein_threshold
            - maximum levenshetein distance the read 
              could have  
    """
    # keep tracks of filtering 
    total_read = -1
    filtered_read = -1
    remained_read = -1
    txt = ("{filt} reads are filtered from {total} reads of {R}, " + 
           "with a levenshtein matching length of {lev_len} and " +
           "threshold of {lev_thres}, " + 
           "resulting in {remain} reads")
    
    # read file
    if file_type in ["fastq.gz", "fastq", "fq"]:
        
        fname = out_name + "." + file_type
        
        if ".gz" in file_type:
            # input file
            file_in = gzip.open(path, "rt")
            # output file
            file_out = gzip.open(out_dir + fname, "wb")
        else:
            # input file
            file_in = open(path, "r")
            #output file
            file_out = open(out_dir + fname, "w")
            
        
        # start reading the reads
        lines = []
        for line in file_in:
            lines.append(line.rstrip())
            
            if len(lines) == 4:
                # reacord the read
                record = process(lines)
                if total_read == -1:
                    total_read = 1
                else:
                    total_read += 1

                read_id = record["read_id"]
                seq = record["seq"]
                sign = record["optional"]
                qual = record["qual"]

                # check the read 
                if read_type == "R1":       
                    checked = check_R1(seq, lev_matching_length,
                                       lev_threshold)
                if read_type == "R2":
                    checked = check_R2(seq, lev_matching_length,
                                       lev_threshold)
                        
                if checked:
                    # pass the check
                    if remained_read == -1:
                        remained_read = 1
                    else:
                        remained_read += 1

                    # write the read to output 
                    if ".gz" in file_type:
                        file_out.write((read_id + "\n" + seq + "\n" + 
                                        sign + "\n" + qual + "\n").encode())
                    else:
                        file_out.write(read_id + "\n" + seq + "\n" + 
                                       sign + "\n" + qual + "\n")

                else:
                    # faile the check
                    if filtered_read == -1:
                        filtered_read = 1
                    else:
                        filtered_read += 1
                
                # empty the line list 
                lines = []
        
        # close the files
        file_in.close()
        file_out.close()
        
        # print output message
        print(txt.format(filt = filtered_read,
                         total = total_read, 
                         R = read_type,
                         lev_len = lev_matching_length,
                         lev_thres = lev_threshold,
                         remain = remained_read))
        print("finished writing filtered reads to " + out_dir + fname + "\n")
            
    else:
        raise ValueError("File type is not fastq.gz, fastq, or fq")

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
    R1_suffix = "_lev" + str(R1_lev_matching_length) + "thres" + str(R1_lev_threshold)
    R1_fname = "R1" + R1_suffix
    filter_read (R1_path, file_type, "R1",
                 R1_fname, out_dir,
             R1_lev_matching_length,
             R1_lev_threshold)

    # process read 2
    R2_suffix = "_lev" + str(R2_lev_matching_length) + "thres" + str(R2_lev_threshold)
    R2_fname = "R2" + R2_suffix
    filter_read (R2_path, file_type, "R2",
                 R2_fname, out_dir,
                 R2_lev_matching_length,
                 R2_lev_threshold)
    
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