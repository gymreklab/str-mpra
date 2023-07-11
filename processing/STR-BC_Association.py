#!/usr/bin/env python3
"""
script for performing STR-barcode association 
v3: directly read in .bam file, filter reads while loading
"""

# Imports 
import argparse
import os
import sys
import copy
import gzip
import pysam
import utils

import numpy as np
import pandas as pd

from cigar import Cigar

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

# Helper functions
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

def getargs():
    parser = argparse.ArgumentParser()    
    # input file & output directory 
    inout_group = parser.add_argument_group("Input/Output")
    inout_group.add_argument("--bam", help="Aligned filtered.bam from previous pre-processing",
                             type=str, required=True)
    inout_group.add_argument("--outdir", help="Path to output directory",
                             type=str, required=True) 
    # filter related value 
    filter_group = parser.add_argument_group("Filter related")
    filter_group.add_argument("--len", help="Expected read length",
                              type=int, required=True)
    filter_group.add_argument("--occurrence", help="Minimum required occurence for a unique STR-BC pair",
                              type=int, default=1)
    filter_group.add_argument("--minBarcode", help="Minimum number of unique barcodes required to be associated per STR",
                              type=int, default=1)
    # get argument
    args = parser.parse_args()
    
    return args

def main(args):    
    # check file existence 
    if not os.path.exists(args.bam):
        print("Error: %s does not exist"%args.bam)
        return 1
    
    # checking if outdir exists, if not, create the outdir
    if not os.path.exists(os.path.dirname(args.outdir)):
         os.mkdir(os.path.dirname(args.outdir))    
    
    # check if summary.csv exists, if not, create a new summary file 
    summ_fname = os.path.join(args.outdir, "summary_STRBC_association.csv")
    if not os.path.exists(summ_fname):
        sum_file = open(summ_fname, "w")
    else:
        sum_file = open(summ_fname, "a")
    
    # Load data
    utils.MSG("Loading raw STR-BC associations")
    bc_str_df = load_bam(args.bam, args.len, sum_file)
    bc_str_df.to_csv(os.path.join(args.outdir, "raw_association.tsv"), \
        sep="\t", index=False)

    # Filter 1: remove BCs corresponding to >1 STR
    # TODO

    # Filter 2: BC occurrence
    # TODO

    # Filter 3: Num. BCs per STR
    # TODO

    # Write output file
    bc_str_df.to_csv(os.path.join(args.outdir, "association.tsv"), \
        sep="\t", index=False)

    sum_file.close()
    utils.MSG("finished writing association.tsv \n")
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