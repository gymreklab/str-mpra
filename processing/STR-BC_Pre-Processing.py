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

def filter_read (path_R1, path_R2,
                 file_type, out_dir, 
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
        4. out_dir
            - output directory 
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
    total_pair = -1
    filtered_pair = -1
    remained_pair = -1
    
    # print message
    txt = ("{filt} pairs of reads are filtered from {total} pairs of reads, " + 
           "with a levenshtein matching length of {R1_lev_len} and " +
           "threshold of {R1_lev_thres} for read 1, " + 
           "and a levenshtein matching length of {R2_lev_len} and " +
           "threshold of {R2_lev_thres} for read 2, " +
           "resulting in {remain} pairs of reads({percent}%).")
    
    # read file
    if file_type not in ["fastq.gz", "fastq", "fq", "fq.gz"]:
        raise ValueError("File type is not fastq.gz, fq.gz, fastq, or fq")
      
    else:        
        fname = os.path.join(out_dir, "filtered" + ".fq.gz")
        
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
                if total_pair == -1:
                    total_pair = 1
                else:
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
                    # pass the check
                    if remained_pair == -1:
                        remained_pair = 1
                    else:
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
                    if filtered_pair == -1:
                        filtered_pair = 1
                    else:
                        filtered_pair += 1
                
                # empty the line list 
                lines = 0
                lines_R1 = []
                lines_R2 = []
        
        # print output message
        percent_remain = "{:.2f}".format((remained_pair/total_pair)*100)
        utils.MSG(txt.format(filt = filtered_pair,
                         total = total_pair, 
                         R1_lev_len = R1_lev_matching_length,
                         R1_lev_thres = R1_lev_threshold,
                         R2_lev_len = R2_lev_matching_length,
                         R2_lev_thres = R2_lev_threshold,
                         remain = remained_pair,
                         percent = percent_remain))
        utils.MSG("finished writing filtered reads to " + out_dir + fname + "\n")
        
        # write in summary 
        sum_file.write("total," + str(total_pair) + "\n")
        sum_file.write("remain," + str(remained_pair) + "\n")
        sum_file.write("% of total," + str(percent_remain) + "\n")
        
        # close the files
        file_R1.close()
        file_R2.close()
        file_out.close()
        sum_file.close()

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
                         /usr/local/bin/samtools view -bS")
                        .format(bwa_ref = ref_path, 
                                reads = read2_path), 
                        shell=True, stdout=bam_file)
        
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
    
    # checking if out_dir exists, if not, create the out_dir
    if not os.path.exists(os.path.dirname(args.outdir)):
         os.mkdir(os.path.dirname(args.outdir))
            
    # create summary.csv file 
    sum_file = open(os.path.join(args.outdir, "summary.csv"), "w")
    
    # process reads 
    utils.MSG("start filtering reads...")
    filter_read(args.read1, args.read2,
                 file_type, args.outdir, 
                 args.R1_match, args.R1_thres, 
                 args.R2_match, args.R2_thres,
                 sum_file)
    utils.MSG("filtering done")

    # alignment on read 2
    out_bam = os.path.join(args.outdir, "filtered.bam")
    utils.MSG("start aligning filtered reads to " + args.bwaref)
    bwamem_alignment(args.bwaref, os.path.join(args.outdir, "filtered." + file_type), out_bam)
    utils.MSG("alignment done")
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