#!/usr/bin/env python3
"""
script for performing STR-barcode association 
v3: directly read in .bam file, filter reads while loading
"""

# Imports 
import os
import sys
import copy
import gzip
import pysam
import argparse
import matplotlib

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from cigar import Cigar

# set plotting style
sns.set(font_scale=2, style="ticks")
plt.rcParams['figure.figsize'] = (60, 12)

# Allow making plots even with no x-forward
matplotlib.use('Agg')

# Allow plots to be editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Helper functions
def cigar_check (cigar_string, R2_length):
    
    """
    check the cigar stirng for qc
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
        check_result = "failed"
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
                    check_result = "failed"
                
        else:
            # cigar string does not start and end with 
            # matching sequence that pass corresponding 
            # threshold, filtered
            check_result = "failed"
    
    return check_result

def num_STR (BC_STR_dict):
    """
    count the number of total STRs in the input barcode-STR dictionary
    """ 
    STRs = set()
    for barcode in BC_STR_dict:
        for STR in BC_STR_dict[barcode]:
            STRs.add(STR)
    
    return len(STRs)

def remove_barcode (BC_STR_dict):
    """
    remove barcodes from the input barcode-STR dictionary that
    are not associate with one and only one STR
    """
    for barcode in BC_STR_dict.copy():
        if len(BC_STR_dict[barcode]) != 1:
            BC_STR_dict.pop(barcode)
            
    return BC_STR_dict
            
def filt_occurrence (BC_STR_dict, thres):
    """
    filter BC-STR pairs from the input barcode-STR dictionary that 
    fail to pass the occurrence threshold
    """
    for barcode in BC_STR_dict.copy():
        for STR in BC_STR_dict[barcode].copy():
            occurrence = BC_STR_dict[barcode][STR]
            if occurrence < thres:
                BC_STR_dict[barcode].pop(STR)
    
    remove_barcode(BC_STR_dict)
    
def load_bam (bam_path, expected_length):
    # output msg
    read_filt_txt = ("Out of {total} reads, {filt} reads that are either " +
                     "less than {expect_len}, and/or does not have " +
                     "a matching STR, and/or fail the cigar check are filtered, " +
                     "resulting in {remain} reads ({percent} %), among which "+
                     "{perfect} reads have a perfect cigar string \n")
    
    # keep track of read filtering 
    num_filt = -1
    remained = -1
    total_read = -1
    num_perfect = -1
    perfect_cig = str(expected_length) + "M"

    # store BC and STR information in the format of 
    # {barcode: {STR: occurrence}} 
    BC_STRs = {}
#     # create a intermediate count.tsv file 
#     count = open(out_dir + "unfiltered_pair_count.tsv", mode="w")
    
    # load aligned reads
    bam_file = pysam.AlignmentFile(bam_path, mode = 'rb')
    bam_iter = bam_file.fetch(until_eof = True)

    for read in bam_iter:

        keep_read = False
        # keep track of total reads
        if total_read == -1:
            total_read = 1
        else:
            total_read += 1

        # obtain read info 
        read_id = str(read.qname)
        barcode = read_id[-20:]
        STR = str(read.reference_name)
        cigar_string = str(read.cigarstring)
        sequence = str(read.query)
        length = len(sequence)

        # keep read that is the expected read length 
        # and does find a matching STR 
        if (length == expected_length) & ("_STR_" in STR):

            # keep read that pass the cigar check
            if cigar_check(cigar_string, expected_length) != "failed":
                keep_read = True

                if cigar_string == perfect_cig:
                    if num_perfect == -1:
                        num_perfect = 1
                    else:
                        num_perfect += 1

        if keep_read:
            # keep track of remaining read
            if remained == -1:
                remained = 1
            else:
                remained += 1

            # count BC and STR
            if barcode not in BC_STRs:
                BC_STRs[barcode] = {STR: 1}
            else:
                if STR not in BC_STRs[barcode]:
                    BC_STRs[barcode][STR] = 1
                else:
                    BC_STRs[barcode][STR] += 1

        else:
            # keep track of filtered read
            if num_filt == -1:
                num_filt = 1
            else:
                num_filt += 1

    bam_file.close()
    print(read_filt_txt.format(total=total_read, filt=num_filt,
                               expect_len=expected_length, 
                               remain=remained,
                               percent=("{:.2f}".format((remained/total_read)*100)),
                               perfect=num_perfect),
          flush=True)
    
#     print("start creating the unfiltered " +
#           "BC-STR pair count matrix...", flush=True)
#     for BC in BC_STRs:
        
#         for STR in BC_STRs[BC]:
#             out = "{pair}\t{occur}\n"
#             count.write(out.format(pair=tuple([BC, STR]),
#                                    occur=BC_STRs[BC][STR]))
        
#     count.close()
#     print("finish creating the unfiltered " +
#           "BC-STR pair count matrix \n", flush=True)
        
    
#     BC_STR_df = pd.DataFrame.from_dict(BC_STRs.items())
#     BC_STR_df.columns = ["BC_STR", "occurrence"]
#     BC_STR_df[['barcode', 'STR']] = pd.DataFrame(BC_STR_df['BC_STR'].tolist(),
#                                                  index=BC_STR_df.index)    
    
    return BC_STRs

def filter_BC_STR (BC_STR_dict, occurrence_thres):
    BC_STRs = copy.deepcopy(BC_STR_dict)
    
    #output msg    
    init_BC_STR_txt = ("{init_barcode} barcodes are captured initially, " +
                       "associating with a total of {init_STR} STRs \n")
    remove_dup_txt = ("{removed} barcodes are found to be associated with, " +
                      "multiple STRs, after removal, {cur_barcode} barcodes " +
                      "are found to be associated with {cur_STR} STRs \n")
    filt_ocur_txt = ("{filtered} unique BC-STR pair is found to occur less than " +
                     "{threshold}, after removal, {fin_barcode} barcodes " +
                      "are found to be associated with {fin_STR} STRs \n")
    
    # record the initial count
    num_init_barcode = len(set(BC_STRs.keys()))
    num_init_STR = num_STR(BC_STRs)
    
    print(init_BC_STR_txt.format(init_barcode=num_init_barcode,
                                 init_STR=num_init_STR),
          flush=True)
    
    # remove barcode associate with multiple STRs
    print("start removing barcodes associate with multiple STRs...",
          flush=True)
    out_dict = remove_barcode(BC_STRs)
    num_cur_barcode = len(set(out_dict.keys()))
    removed_barcode = num_init_barcode - num_cur_barcode
    num_cur_STR = num_STR(out_dict)
    
    print(remove_dup_txt.format(removed=removed_barcode, 
                                cur_barcode=num_cur_barcode,
                                cur_STR=num_cur_STR))

    # filter on BC-STR pair occurrence
    print("start filtering on BC-STR pair occurrence...",
          flush=True)
    filt_occurrence(out_dict, occurrence_thres)
    num_fin_barcode = len(set(out_dict.keys()))
    removed_barcode = num_cur_barcode - num_fin_barcode
    num_fin_STR = num_STR(out_dict)
    
    print(filt_ocur_txt.format(filtered=removed_barcode,
                               threshold=occurrence_thres,
                               fin_barcode=num_fin_barcode,
                               fin_STR=num_fin_STR))
    
    return out_dict

def countplot_STRBC_occurrence (BC_STR_dict, out_dir, fig_suffix=None):
    plt.figure()
    # generate occurrence list
    occurrences = []
    for barcode in BC_STR_dict:
        for STR in BC_STR_dict[barcode]:
            occurrences.append(BC_STR_dict[barcode][STR])

    oc = sns.countplot(x=occurrences);
    plt.title("How many BC-STR pairs occurred multiple time?")
    plt.xticks(fontsize=8)

#     oc2 = plt.axes([0.3, 0.3, 0.55, 0.55])
#     sns.countplot(x=occurrences, ax = oc2);
#     oc2.set_title('zoom')
#     oc2.set_xlabel(None)
#     oc2.set_ylabel(None)
#     oc2.set_ylim([0,120]);

    if fig_suffix is None:
        plt.savefig(out_dir + 'BC-STR_occurrence.png')
    else:
        plt.savefig(out_dir + fig_suffix + 'BC-STR_occurrence.png')
            
def STR_count (BC_STR_dict, out_dir=False):
    """
    create a dictionary of {STR: num of barcode associated}
    """
    
    STR_count = {}
    
    for barcode in BC_STR_dict:
        for STR in BC_STR_dict[barcode]:
            if STR not in STR_count:
                STR_count[STR] = 1
            else:
                STR_count[STR] += 1
    
    if out_dir:
        print("start writing the number of barcode associated with " +
              "each STR to BC_per_STR.tsv in the format of:\n " + 
              "STR\t num barcode associated with this STR ",
              flush=True)

        path = out_dir + "BC_per_STR.tsv"
        file = open(path, "w")

        for STR in STR_count:
            file.write(STR + "\t" + str(STR_count[STR]) + "\n")

        file.close()
        print("finished writing BC_per_STR.tsv\n")
                
    return STR_count

def filt_num_barcode (BC_STR_dict, STR_count_dict, barcode_thres):
    BC_STRs = copy.deepcopy(BC_STR_dict)
    
    """
    filter BC-STR pairs from the input barcode-STR dictionary that 
    fail to pass the occurrence threshold
    """
    
    # record the initial count
    num_init_barcode = len(set(BC_STRs.keys()))
    num_init_STR = num_STR(BC_STRs)
    
    # filtering
    for barcode in BC_STRs.copy():
        for STR in BC_STRs[barcode].copy():
            num_barcode = STR_count_dict[STR]
            if num_barcode < barcode_thres:
                BC_STRs[barcode].pop(STR)
    
    out_dict = remove_barcode(BC_STRs)
    
    #output msg
    filt_bc_txt = ("{filtered} STR is found to have less than " +
                   "{threshold} barcodes associated with, after removal, "+ 
                   "{fin_barcode} barcodes are found to be associated with {fin_STR} STRs \n")
    # filter on BC-STR pair occurrence
    print("start filtering on num barcode associated per STR...",
          flush=True)
    num_fin_barcode = len(set(out_dict.keys()))
    removed_barcode = num_init_barcode - num_fin_barcode
    num_fin_STR = num_STR(out_dict)
    
    print(filt_bc_txt.format(filtered=removed_barcode,
                             threshold=barcode_thres,
                             fin_barcode=num_fin_barcode,
                             fin_STR=num_fin_STR))
    
    return out_dict

def countplot_multiBC (STR_count_dict, out_dir, fig_suffix=None):
    
    plt.figure()
    STR_df = pd.DataFrame.from_dict(data=STR_count_dict.items())
    STR_df.columns=["STR", "BC_count"]

    # output msg
    txt = ("The maximum barcode associate with a unique STR is {maxSTR}, " +
           "and there are {maxSTR_count} such unique STR, " +
           "the minimum barcode associate with a unique STR is {minSTR}, " +
           "and there are {minSTR_count} such unique STR")

    # maximum barcode associated with 1 unique STR
    STR_maxBC = STR_df.BC_count.max()
    # number of such STR
    STR_maxBC_count = STR_df["BC_count"].value_counts()[STR_maxBC]

    # minimum barcode associated with 1 unique STR
    STR_minBC = STR_df.BC_count.min()
    # number of such STR
    STR_minBC_count = STR_df["BC_count"].value_counts()[STR_minBC]

    # plotting 
    mbc = sns.countplot(data=STR_df, x="BC_count")
    mbc.set(title="Number of unique barcodes associated per unique STR");
    mbc.set_xlabel("Number of unique barcode")
    plt.xticks(fontsize=8)

#     mbc2 = plt.axes([0.3, 0.3, 0.55, 0.55])
#     sns.countplot(data=STR_df, x="BC_count", ax = mbc2)
#     mbc2.set_title('zoom')
#     mbc2.set_xlabel(None)
#     mbc2.set_ylabel(None)
#     mbc2.set_ylim([0,120]);
    
    if fig_suffix is None:
        plt.savefig(out_dir + 'BC_per_STR.png')
    else:
        plt.savefig(out_dir + fig_suffix + 'BC_per_STR.png')

    print(txt.format(maxSTR = STR_maxBC,
                     maxSTR_count = STR_maxBC_count,
                     minSTR = STR_minBC,
                     minSTR_count = STR_minBC_count), 
          flush=True)
        
def output_count_and_association (BC_STR_dict, out_dir):
    
    print("start writing to association.tsv in the format of:\n " +
          "barcode\t STR\t occurrence",
          flush=True)
    
    path = out_dir + "association.tsv"
    file = open(path, "w")
    
    for barcode in BC_STR_dict:
        for STR in BC_STR_dict[barcode]:
            occurrence = BC_STR_dict[barcode][STR]
            out = "{barcode}\t{STR}\t{occur}\n"
            file.write(out.format(barcode=barcode,
                                  STR=STR,
                                  occur=occurrence))   

    file.close()
    print("finished writing association.tsv \n")

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
                              type=int, required=True)
    filter_group.add_argument("--minBarcode", help="Minimum number of unique barcodes required to be associated per STR",
                              type=int, required=True)
    
    # optional suffix for plot
    opt_group = parser.add_argument_group("Optional suffix for plots")
    opt_group.add_argument("--occurCount", help="Suffix for STR-BC pair occurrence count plot",
                           type=str, default=None)
    opt_group.add_argument("--multiBC", help="Suffix for count plot of number of barcode associate with one unique STR",
                           type=str, default=None)
    
    # get argument
    args = parser.parse_args()
    
    return args

def main(args):
    
    # required input and output path
    bam_path = args.bam
    out_dir = args.outdir
    
    # required filter related parameter 
    expected_length = args.len
    occurrence_thres = args.occurrence 
    barcode_thres = args.minBarcode

    # check file existence 
    if not os.path.exists(bam_path):
        common.WARNING("Error: %s does not exist"%bam_path)
        return 1
    
    # checking if out_dir exists, if not, create the out_dir
    if not os.path.exists(os.path.dirname(out_dir)):
         os.makedir(os.path.dirname(out_dir))    
    
    # optional suffix for plots
    STRBC_occurrence_plot_suffix = args.occurCount
    multiBC_plot_suffix = args.multiBC
    
    # load data 
    BC_STR = load_bam(bam_path, expected_length)

    # filtering
    filt_occurrence = filter_BC_STR(BC_STR, occurrence_thres)
    init_STR = STR_count(filt_occurrence)
    fin_BC_STR = filt_num_barcode(filt_occurrence, init_STR, barcode_thres)
    
    # output
    fin_STR = STR_count(fin_BC_STR, out_dir)
    output_count_and_association (fin_BC_STR, out_dir)
    
    # plot occurrence distribution 
    print("start plotting occurrence distribution...",
          flush=True)
    countplot_STRBC_occurrence (fin_BC_STR, out_dir, STRBC_occurrence_plot_suffix)
    print("finished plotting occurrence distribution\n",
          flush=True)
    
    # plot num barcode per STR
    print("start plotting the number of barcode that is associated with " +
          "a unique STR...",
          flush=True)
    countplot_multiBC (fin_STR, out_dir, multiBC_plot_suffix)
    print("finished plotting number of barcode per STR\n",
          flush=True)

    print("BC-STR association done")
    
    
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