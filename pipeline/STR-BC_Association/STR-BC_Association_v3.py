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
plt.rcParams['figure.figsize'] = (20, 12)

# Allow making plots even with no x-forward
matplotlib.use('Agg')

# Allow plots to be editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Helper functions
def load_R1 (R1_path):
    
    # load filtered R1 
    reads = {"read_id":[],
             "sequence":[],
             "quality":[]}
    
    if ".gz" in R1_path:
        with gzip.open(R1_path, "rt") as file:
            df = pd.read_csv(file, sep='\n', header=None)
            filt_R1 = pd.DataFrame(df.values.reshape(-1, 4))
            filt_R1.columns=['read_id', "R1", '+', 'qual']
            
        file.close()

    else:
            df = pd.read_csv(R1_path, sep='\n', header=None)
            filt_R1 = pd.DataFrame(df.values.reshape(-1, 4))
            filt_R1.columns=['read_id', "R1", '+', 'qual']
    
    # only keep read id and sequence
    filt_R1 = filt_R1[["read_id", "R1"]]
    
    # get seq length
    filt_R1["length"] = filt_R1.R1.str.len()

    # modify read_id for easier merging 
    filt_R1["read_id"] = filt_R1.read_id.str[1:]
    filt_R1_id = filt_R1["read_id"].str.split(" ", expand=True)
    filt_R1["read_id"] = filt_R1_id[0]

    # obtain barcode for each read
    filt_R1["barcode"] = filt_R1.R1.str[:20]
    
    return filt_R1

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

def load_R2 (R2_path, R2_length):
    
    num_filt = 0
    total_read = 0
    num_perfect = 0
    perfect_cig = str(R2_length) + "M"
    
    # load aligned R2
    bam_file = pysam.AlignmentFile(R2_path, mode = 'rb')
    bam_iter = bam_file.fetch(until_eof = True)

    bam_dict = {
        "read_id": [],
        "STR": [],
        "CIGAR":[],
        "R2":[],
        "length":[]
    }
    for read in bam_iter:
        
        keep_read = False
        total_read += 1
        
        read_id = str(read.qname)
        STR = str(read.reference_name)
        cigar_string = str(read.cigarstring)
        sequence = str(read.query)
        length = len(sequence)
        
        # keep read that is the expected read length 
        # and does find a matching STR 
        if (length == R2_length) & ("_STR_" in STR):
            
            # keep read that pass the cigar check
            if cigar_check(cigar_string, R2_length) != "failed":
                keep_read = True
                
                if cigar_string == perfect_cig:
                    num_perfect += 1
        
        if keep_read:       
            bam_dict["read_id"].append(read_id)
            bam_dict["STR"].append(STR)
            bam_dict["CIGAR"].append(cigar_string)
            bam_dict["R2"].append(sequence)
            bam_dict["length"].append(length)
            
        else:
            num_filt += 1
        
    bam_file.close()
    
    aln_R2 = pd.DataFrame(bam_dict)
    
    txt = ("Out of {total} read 2, {filt} reads that are either " +
           "less than {expect_len}, and/or does not have " +
           "a matching STR, and/or fail the cigar check are filtered, " +
           "resulting in {remain} reads({percent}%), among which "+
           "{perfect} reads have a perfect cigar string \n")
    
    print(txt.format(total=total_read, filt=num_filt,
                     expect_len=R2_length, 
                     remain=len(aln_R2),
                     percent=("{:.2f}".format((len(aln_R2)/total_read)*100)),
                     perfect=num_perfect), 
          flush=True)
    
    return aln_R2

def BC_check (in_df):

    STR_BCs = []
    BCs = []
    STRs = []
    occurrence = []

    for STR, barcode in zip(in_df.STR, in_df.barcode):
        STR_BC = tuple([STR, barcode])

        if STR_BC not in STR_BCs:
            STR_BCs.append(STR_BC)
            BCs.append(barcode)
            STRs.append(STR)
            occurrence.append(1)
        else:
            occurrence[STR_BCs.index(STR_BC)] += 1

    check = pd.DataFrame({"STR_BC": STR_BCs,
                          "barcode": BCs, 
                          "STR": STRs,
                          "occurrence": occurrence})
    
    return check

def remove_dup_BC (check):
    duplicate_barcode = check[check.duplicated('barcode')]
    remove_dup = set(duplicate_barcode.barcode)
    print(str(len(remove_dup)) + " unique barcodes are found to be " +
          "associated with multiple STRs", 
          flush=True)
    print("remove such barcode \n", 
          flush=True)
    
    association_df = check[~check["barcode"].isin(remove_dup)]

    return association_df

def filt_occurrence (association_df, thres):
    before = len(association_df)
    df = copy.deepcopy(association_df)
    
    df = df[df["occurrence"] >= thres]
    print(("{filtered} unique STR-BC pair is found to occur less than {threshold}")
          .format(filtered=(before - len(df)),
                  threshold=thres), 
          flush=True)
    print("remove such STR-BC pair \n",
          flush=True)
    
    return df

def countplot_STRBC_occurrence (association_df, out_dir, fig_suffix=None):
    plt.figure()
    oc = sns.countplot(data=association_df, x="occurrence");
    plt.title("How many STR-BC pairs occurred multiple time?")
    
    oc2 = plt.axes([0.3, 0.3, 0.55, 0.55])
    sns.countplot(data=association_df, x="occurrence", ax = oc2)
    oc2.set_title('zoom')
    oc2.set_xlabel(None)
    oc2.set_ylabel(None)
    oc2.set_ylim([0,120]);
    
    if fig_suffix is None:
        plt.savefig(out_dir + 'STR-BC_occurrence.png')
    else:
        plt.savefig(out_dir + fig_suffix + 'STR-BC_occurrence.png')
        
def countplot_multiBC (association_df, out_dir, fig_suffix=None):
    plt.figure()
    STR_count = association_df.STR.value_counts().to_frame().reset_index()
    STR_count.columns = ["STR", "count"]
    STR_count["count"].value_counts()
    # output information
    txt = ("The maximum barcode associate with a unique STR is {maxSTR}, " +
           "and there are {maxSTR_count} such unique STR, " +
           "the minimum barcode associate with a unique STR is {minSTR}, " +
           "and there are {minSTR_count} such unique STR \n")
    
    # maximum barcode associated with 1 unique STR
    STR_maxBC = association_df.STR.value_counts().max()
    # number of such STR
    STR_maxBC_count = STR_count["count"].value_counts()[STR_maxBC]
    
    # minimum barcode associated with 1 unique STR
    STR_minBC = association_df.STR.value_counts().min()
    # number of such STR
    STR_minBC_count = STR_count["count"].value_counts()[STR_minBC]
    
    # plotting 
    mbc = sns.countplot(data=STR_count, x="count")
    mbc.set(title="After filtering STR-BC pairs with low occurrence,\nhow many barcode is associated with a unique STR?");
    mbc.set_xlabel("Number of barcode associate with unique STR")
    
    mbc2 = plt.axes([0.3, 0.3, 0.55, 0.55])
    sns.countplot(data=STR_count, x="count", ax = mbc2)
    mbc2.set_title('zoom')
    mbc2.set_xlabel(None)
    mbc2.set_ylabel(None)
    mbc2.set_ylim([0,120]);
    
    if fig_suffix is None:
        plt.savefig(out_dir + 'multi-BC_count.png')
    else:
        plt.savefig(out_dir + fig_suffix + 'multi-BC_count.png')

    print(txt.format(maxSTR = STR_maxBC,
                     maxSTR_count = STR_maxBC_count,
                     minSTR = STR_minBC,
                     minSTR_count = STR_minBC_count), 
          flush=True)
        
def output_association (association_df, out_dir):

    num_STR = len(association_df.STR.unique())
    num_BC = len(association_df.barcode.unique())
                  
    path = out_dir + "association.tsv"
    file = open(path, "w")
    print("start writing association.tsv", 
          flush=True)
    
    for STR_BC in association_df.STR_BC:
        barcode = STR_BC[1]
        STR = STR_BC[0]
        
        file.write(barcode + "\t" + STR + "\n")   

    file.close()
        
    print(str(num_BC) + " unique barcode is matched to " + 
          str(num_STR) + " unique STR \n", 
          flush=True)

def getargs():
    parser = argparse.ArgumentParser()
    
    # input file & output directory 
    inout_group = parser.add_argument_group("Input/Output")
    inout_group.add_argument("--filtR1", help="Filtered read1 from BC_read_processing.py",
                             type=str, required=True)
    inout_group.add_argument("--alnR2", help="Aligned read2.bam from BC_read_processing.py",
                             type=str, required=True)
    inout_group.add_argument("--outdir", help="Path to output directory",
                             type=str, required=True)
    
    # filter related value 
    filter_group = parser.add_argument_group("Filter related")
    filter_group.add_argument("--lenR1", help="Expected read length of read1",
                              type=int, required=True)
    filter_group.add_argument("--lenR2", help="Expected read length of read2",
                              type=int, required=True)
    filter_group.add_argument("--occurrence", help="Minimum required occurence for a unique STR-BC pair",
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
    filt_R1_path = args.filtR1
    aln_R2_path = args.alnR2
    out_dir = args.outdir
    
    # required filter related parameter 
    R1_length = args.lenR1
    R2_length = args.lenR2
    occurrence_thres = args.occurrence 

    # check file existence 
    if not os.path.exists(filt_R1_path):
        print("Error: %s does not exist"%filt_R1_path)
        return 1
    
    if not os.path.exists(aln_R2_path):
        print("Error: %s does not exist"%aln_R2_path)
        return 1
    
    if not os.path.exists(os.path.dirname(os.path.abspath(out_dir))):
        common.WARNING("Error: The output directory {outdir} does not exist"
                       .format(outdir=out_dir))
        return 1    
    
    # optional suffix for plots
    STRBC_occurrence_plot_suffix = args.occurCount
    multiBC_plot_suffix = args.multiBC
    
    # load data 
    # load R1
    filt_R1 = load_R1(filt_R1_path)
    # load R2
    aln_R2 = load_R2(aln_R2_path, R2_length)

    # merge reamining read 2(STR) with read 1(BC) by read ID on read 2
    merge = pd.merge(left=aln_R2, right=filt_R1,
                     how='left', left_on="read_id", right_on="read_id")
    # the amount of R2 without a matching R1
    R2_count = len(merge)
    R2_no_R1 = int(merge.isnull().sum().R1)
    print("Out of " + str(R2_count) + " read2s, " + str(R2_no_R1) + 
          "(" + str("{:.2f}".format((R2_no_R1/R2_count)*100)) + "%)" +
          " of the reads does not have a read1 with the same read_id and are removed\n")
    # drop the null value
    merge = merge.dropna()

    # filter BC associate with multiple STRs
    check = BC_check(merge)
    association_df = remove_dup_BC(check)

    # plot occurrence distribution 
    countplot_STRBC_occurrence(association_df, out_dir, STRBC_occurrence_plot_suffix)

    # filter STR-BC pair with occurrence less than occurrence_thres
    filt_occurrence(association_df, occurrence_thres)

    # plot how many STR is associated with multiple BC
    countplot_multiBC(association_df, out_dir, multiBC_plot_suffix)
    
    # output
    output_association(association_df, out_dir)
    
    
    
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