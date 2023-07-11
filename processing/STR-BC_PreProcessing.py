#!/usr/bin/env python3
"""
script for pre-processing of STR-BC assocition 
"""

# imports
import argparse
import copy
import gzip
import Levenshtein
import os
import numpy as np
import pandas as pd
import pysam
import subprocess
import sys
import utils

# helper functions
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

def cigar_check (cigar_string, R2_length):
    
    """
    check the cigar string for qc
    parameter: 
      cigar_string - the cigar string from col 6 of sam file 
    output:
      if cigar string pass the check, 
      return check_result = cigar_string
      if not, return check_result = "failed"
    """
    
    # internal thresholds used
    expected_length = 4 
    start_M_thres = 25 # min num bp match at start
    end_M_thres = 20 # min num bp match at end 
    indel_thres = 2 # max num bp of indel allowance 
    
    # check
    perfect_match = str(R2_length) + "M"
    if cigar_string == perfect_match:
        # no need to check if cigar string 
        # is perfetc match
        return cigar_string
    
    parse_cigar = list(Cigar(cigar_string).items())
    
    if len(parse_cigar) > expected_length:
        # long cigar string, filtered 
        return "failed"
    else:
        start = parse_cigar[0]
        end = parse_cigar[-1]
        mid = parse_cigar[1:-1]
        
        if ((start[1]=="M" and start[0] >= start_M_thres) and
            (end[1]=="M" and end[0] >= end_M_thres)):
            
            for sub_cigar in mid:
                if ((sub_cigar[1] in ["I", "D"] ) and
                    (sub_cigar[0] < indel_thres)):
                    check_result = cigar_string
                else:
                    # cigar string either have modification 
                    # other than indel, or the indel is larger
                    # than the corresponding threshold
                    return "failed"    
        else:
            # cigar string does not start and end with 
            # matching sequence that pass corresponding 
            # threshold, filtered
            return "failed"
    
    return check_result

def load_bam(bam_path, expected_length, sum_file):
    data_barcodes = []
    data_strs = []

    # keep track of read filtering 
    num_filt = 0
    remained = 0
    total_read = 0
    num_perfect = 0
    perfect_cig = str(expected_length) + "M"
    
    # load aligned reads
    bam_file = pysam.AlignmentFile(bam_path, mode = 'rb')
    bam_iter = bam_file.fetch(until_eof = True)

    for read in bam_iter:
        keep_read = False
        total_read += 1

        # obtain read info 
        barcode = str(read.qname)[-20:]
        STR = str(read.reference_name)
        cigar_string = str(read.cigarstring)
        length = len(read.query)

        # keep read that is the expected read length 
        # and does find a matching STR 
        if (length == expected_length) & ("STR" in STR):
            if cigar_check(cigar_string, expected_length) != "failed":
                keep_read = True
                if cigar_string == perfect_cig: num_perfect += 1
        if keep_read:
            remained += 1
            # count BC and STR
            data_barcodes.append(barcode)
            data_strs.append(STR)
        else:
            num_filt += 1

    bam_file.close()
    percent_remain = "{:.2f}".format((remained/total_read)*100)
    # write in summary 
    sum_file.write("aligned and pass association qc," + str(remained) + "\n")
    sum_file.write("% of total," + str(percent_remain) + "\n")
    sum_file.write("perfect cigar," + str(num_perfect) + "\n")
    df =pd.DataFrame({"barcode": data_barcodes,
        "STR": data_strs})
    df["count"] = 1
    df = df.groupby(["barcode","STR"], as_index=False).agg({"count": np.sum})
    return df

def filter_read (path_R1, path_R2,
                 file_type, out_prefix, 
                 R1_lev_matching_length,
                 R1_lev_threshold, 
                 R2_lev_matching_length,
                 R2_lev_threshold,
                 sum_file):
    """
    to filter reads and write the filter read 
    to a new read file of the same file type
    parameter:
        1. path_R#
            - path to read file
        2. file_type
            - type of read file, fastq.gz or fastq or fq
        3. read_type
            - type of read, R1 or R2
        4. out_prefix
            - prefix to name output files
        5. R#_levenshetein_matching_length
            - the length of the sequence that will be 
              used to calculate levenshetein distancing 
        6. R#_levenshetein_threshold
            - maximum levenshetein distance the read 
              could have  
        7. sum_file 
            - summary statistic file 
    """
    
    # keep tracks of filtering 
    total_pair = 0
    filtered_pair = 0
    remained_pair = 0
    
    # read file
    if file_type not in ["fastq.gz", "fastq", "fq", "fq.gz"]:
        raise ValueError("File type is not fastq.gz, fq.gz, fastq, or fq")
      
    else:        
        fname = out_prefix + ".filtered.fq.gz"
        
        if ".gz" in file_type:
            # input file
            file_R1 = gzip.open(path_R1, "rt")
            file_R2 = gzip.open(path_R2, "rt")
        else:
            # input file
            file_R1 = open(path_R1, "r")
            file_R2 = open(path_R2, "r")
        file_out = open(fname, "wb")
            
        # start reading the reads
        lines = 0
        lines_R1 = []
        lines_R2 = []
        for line_R1, line_R2 in zip(file_R1, file_R2):
            lines_R1.append(line_R1.rstrip())
            lines_R2.append(line_R2.rstrip())
            lines += 1
            
            if lines == 4:
                
                # record the read
                record_R1 = process(lines_R1)
                record_R2 = process(lines_R2)
                total_pair += 1

                # record R1
                R1_id = record_R1["read_id"]
                R1_seq = record_R1["seq"]
                R1_sign = record_R1["optional"]
                R1_qual = record_R1["qual"]
                
                # record R2
                R2_id = record_R2["read_id"]
                R2_seq = record_R2["seq"]
                R2_sign = record_R2["optional"]
                R2_qual = record_R2["qual"]
                
                # check the read 
                R1_pass = check_R1(R1_seq, R1_lev_matching_length,
                                    R1_lev_threshold)
                R2_pass = check_R2(R2_seq, R2_lev_matching_length,
                                    R2_lev_threshold)
                
                if R1_pass & R2_pass:
                    remained_pair += 1
                    # write the read to output 
                    barcode = R1_seq[0:20]
                    R2_id_lists = list(R2_id.split(" "))
                    R2_id = (R2_id_lists[0] + ":" + barcode +
                             " " + R2_id_lists[1])                    
                    file_out.write((R2_id + "\n" + R2_seq + "\n" + 
                                        R2_sign + "\n" + R2_qual + "\n").encode())
                else:
                    # failed the check
                    filtered_pair += 1
                
                # empty the line list 
                lines = 0
                lines_R1 = []
                lines_R2 = []
        
        percent_remain = "{:.2f}".format((remained_pair/total_pair)*100)        
        # write in summary 
        sum_file.write("Total reads," + str(total_pair) + "\n")
        sum_file.write("Reads remaining," + str(remained_pair) + "\n")
        sum_file.write("Percent of total," + str(percent_remain) + "\n")
        
        # close the files
        file_R1.close()
        file_R2.close()
        file_out.close()

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
        subprocess.call(("bwa mem -t 12 -L 100 -k 8 -O 5 {bwa_ref} {reads} |\
                         /usr/local/bin/samtools view -bS |\
                         samtools sort -O bam")
                        .format(bwa_ref = ref_path, 
                                reads = read2_path), 
                        shell=True, stdout=bam_file, 
                        stderr=open(out_path + ".bwa.log", "w"))
    os.system("samtools index %s"%out_path)

def getargs():
    parser = argparse.ArgumentParser()
    # input file & output directory 
    inout_group = parser.add_argument_group("Input/Output")
    inout_group.add_argument("--read1", help="Path to read1 file", type=str,
                             required=True)
    inout_group.add_argument("--read2", help="Path to read2 file", type=str,
                             required=True)
    inout_group.add_argument("--bwaref", help="Path to ref.fa", type=str,
                             required=True)
    inout_group.add_argument("--outprefix", help="Prefix to name output files", type=str,
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
    filter_group.add_argument("--len", help="Expected read length",
                              type=int, required=True)
    # get argument
    args = parser.parse_args()
    
    return args
    
def InferFileType(fqfile):
    if fqfile.endswith(".fastq.gz"):
        return "fastq.gz"
    elif fqfile.endswith(".fastq"):
        return "fastq"
    elif fqfile.endswith(".fq"):
        return "fq"
    elif fqfile.endswith(".fq.gz"):
        return "fq.gz"

def main(args):
    # input/output
    file_type = InferFileType(args.read1)

    # check file existence 
    if not os.path.exists(args.read1):
        utils.MSG("Error: %s does not exist"%args.read1)
        return 1
    
    if not os.path.exists(args.read2):
        utils.MSG("Error: %s does not exist"%args.read2)
        return 1
    
    if not os.path.exists(args.bwaref):
        utils.MSG("Error: %s does not exist"%args.bwaref)
        return 1
    
    # checking if output directory exists, if not, create it
    if not os.path.exists(os.path.dirname(args.outprefix)):
         os.mkdir(os.path.dirname(args.outprefix))
            
    # create summary.csv file 
    sum_file = open(args.outprefix + ".summary.csv", "w")
    
    # Filter reads to get new fq for read2
    utils.MSG("Filter STR-BC Reads...")
    filter_read(args.read1, args.read2,
                 file_type, args.outprefix, 
                 args.R1_match, args.R1_thres, 
                 args.R2_match, args.R2_thres,
                 sum_file)

    # Perform alignment on read 2 only
    utils.MSG("Align filtered read 2 to " + args.bwaref)
    out_bam = args.outprefix + ".filtered.bam"
    bwamem_alignment(args.bwaref, args.outprefix + ".filtered.fq.gz", out_bam)

    # Load BAM to get STR-BC table
    utils.MSG("Loading STR-BC info")
    bc_str_df = load_bam(out_bam, args.len, sum_file)
    bc_str_df.to_csv(args.outprefix +  ".raw_association.tsv", \
        sep="\t", index=False)

    # Filter 1: remove BCs corresponding to >1 STR
    # TODO

    # Filter 2: BC occurrence
    # TODO

    # Filter 3: Num. BCs per STR
    # TODO

    # Write output file
    utils.MSG("Writing final STR-BC association file")
    bc_str_df.to_csv(args.outprefix + ".association.tsv", \
        sep="\t", index=False)
    sum_file.close()
    utils.MSG("Done!")
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