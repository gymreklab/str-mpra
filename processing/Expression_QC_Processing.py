#!/usr/bin/env python3
"""
script for quality checking and processing the expression reads 
which will output raw count matricies, count matricies that 
filtered out those with low read counts, aggregated count matrix, 
scatter plot within and between replicates, and distribution plots. 

can also output filtered.fastq.gz and failed.fastq.gz upon request.
"""

#set up
import os
import sys 
import copy 
import gzip
import argparse
import Levenshtein 
import numpy as np
import pandas as pd 
import scipy.stats as stats

import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
import matplotlib.colors as mcolor

# set plotting style
sns.set(font_scale=2, style="ticks")
plt.rcParams['figure.figsize'] = (40, 12)

# Allow making plots even with no x-forward
matplotlib.use('Agg')

# Allow plots to be editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# helper functions
def process (lines=None):
    
    """
    helper function to process the fastq read 
    """
    
    ks = ['read_id', 'seq', 'optional', 'qual']
    return {k: v for k, v in zip(ks, lines)}

def check_read (sequence, read_len,
                levenshtein_matching_length,
                levenshtein_threshold):
    """
    to check the expression reads
    parameter:
        1. sequence
            - the read sequence
        2. read_len
            - the length of the read
        3. levenshetein_matching_length
            - the length of the sequence that will be 
              used to calculate levenshetein distancing 
        4. levenshetein_threshold
            - maximum levenshetein distance the read 
              could have 
    output:
        return true if the sequence pass the check, 
        return false other wise
    """
    
    # the full identical GFP sequence that could be 
    # used for qc
    qc_seq = "TCTAGAGGTTCGTCGACGCGATTATTATCATTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGG"
    
    # the sqeuence required to be exact match 
    exact_match = qc_seq[:levenshtein_matching_length]
    
    # the sequence being checked 
    # if the input read length is longer than the NovaSeq read length, 
    # capped at the NovaSeq read length
    # otherwise use the shorter read length 
    if read_len >= 101:
        read_len = 101 
        
    seq_check = sequence[20:read_len]
    
    # check for required bp of exact match after barcode
    if seq_check[:levenshtein_matching_length] == exact_match:
        
        # calculate levenshtein distancing of overall sequence
        lev_score = Levenshtein.distance(qc_seq[:read_len-20], 
                                         seq_check)
        # check for leveshtein distancing
        if lev_score > levenshtein_threshold:
            # have more mismtach than allowed 
            # fail the check 
            return False
        
    else:
    # does not have the required bp of exact match after barcode
    # fail the check
        return False
    
    # pass the check
    return True

def filter_read (path, file_type,
                 read_type, read_len, 
                 group_num, out_dir,
                 lev_matching_length,
                 lev_threshold, 
                 sum_dict, export_file=None):
    """
    to filter reads and write the filter read 
    to a new read file of the same file type,
    and generate a barcode count dict during the
    filtering process 
    parameter:
        1. path 
            - path to read file
        2. file_type
            - type of read file, fastq.gz or fastq or fq
        3. read_type
            - type of read, cDNA or gDNA
        4. group_num
            - replicate number of the file 
        5. out_dir
            - output directory 
        6. levenshetein_matching_length
            - the length of the sequence that will be 
              used to calculate levenshetein distancing 
        7. levenshetein_threshold
            - maximum levenshetein distance the read 
              could have  
        8. sum_dict
            - summary dictionary used to record stats 
        9. export_file
            - if not None, output filtered and failed .fastq.gz 
    output:
        barcode_dict
        - a barcode count dict
    """
    # keep tracks of filtering 
    total_read = -1
    filtered_read = -1
    remained_read = -1
    txt = ("{filt} reads are filtered from {total} reads of {R} " + 
           "with a levenshtein matching length of {lev_len} and " +
           "threshold of {lev_thres}, " + 
           "resulting in {remain} reads ({percent}% remains)")
    barcode_dict = {}
    
    # read file
    if file_type in ["fastq.gz", "fastq", "fq"]:
        
        # file naming 
        suffix = "_lev" + str(lev_matching_length) + "thres" + str(lev_threshold)
        out_name = read_type + str(group_num) + suffix
        fname = out_name + "." + file_type
        fail = "failed_" + out_name + "." + file_type
        
        if ".gz" in file_type:
            # input file
            file_in = gzip.open(path, "rt")
            if export_file:
                # output file
                file_out = gzip.open(out_dir + fname, "wb")
                file_fail = gzip.open(out_dir + fail, "wb")
        else:
            # input file
            file_in = open(path, "r")
            if export_file:
                #output file
                file_out = open(out_dir + fname, "w")
                file_fail = open(out_dir + fail, "w")
        
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
                checked = check_read(seq, read_len, 
                                     lev_matching_length,
                                     lev_threshold)
                        
                if checked:
                    # pass the check
                    if remained_read == -1:
                        remained_read = 1
                    else:
                        remained_read += 1
                    
                    # record the barcode 
                    barcode = seq[:20]
                    if barcode in barcode_dict:
                        barcode_dict[barcode] += 1
                    else:
                        barcode_dict[barcode] = 1
                    
                    if export_file:
                        # write the read to output 
                        if ".gz" in file_type:
                            file_out.write((read_id + "\n" + seq + "\n" + 
                                            sign + "\n" + qual + "\n").encode())
                        else:
                            file_out.write(read_id + "\n" + seq + "\n" + 
                                           sign + "\n" + qual + "\n")

                else:
                    # fail the check
                    if filtered_read == -1:
                        filtered_read = 1
                    else:
                        filtered_read += 1
                    
                    if export_file:
                        # output read that fail the qc
                        if ".gz" in file_type:
                            file_fail.write((read_id + "\n" + seq + "\n" + 
                                             sign + "\n" + qual + "\n").encode())
                        else:
                            file_fail.write(read_id + "\n" + seq + "\n" + 
                                            sign + "\n" + qual + "\n")
                
                # empty the line list 
                lines = []
        
        # close the input and output files
        file_in.close()
        if export_file:
            file_out.close()
            file_fail.close()
            print("finished writing filtered reads to " + out_dir + fname + "\n",
                  flush=True)
        
        # print output message
        percent_remained = "{:.2f}".format((remained_read/total_read)*100)
        print(txt.format(filt = filtered_read,
                         total = total_read, 
                         R = read_type + str(group_num),
                         lev_len = lev_matching_length,
                         lev_thres = lev_threshold,
                         remain = remained_read,
                         percent = percent_remained),
              flush=True)
        
        # record the statistic
        sum_dict[read_type + str(group_num)] = {
            "total": total_read,
            "passed_qc": remained_read,
            "percent_passed": percent_remained, 
        }
        
        return barcode_dict
            
    else:
        raise ValueError("File type is not fastq.gz, fastq, or fq")
        
def count_filter (in_df, threshold):
    df = copy.deepcopy(in_df)
    
    for col in df.columns:
        if "gDNA" in col:
            gDNA = col
        if "cDNA" in col:
            cDNA = col
    
    # filter barcode with corressponding count
    # less than threshold 
    df = df[df[gDNA] >= threshold]
    df = df[df[cDNA] >= threshold]
    
    return df

def BC_STR_association (count_filtered_df, out_dir,
                        association, min_barcode,
                        group_num, sum_dict):

    out_msg = ("Out of {total} unique barcode, {no_STR}({percent}%) " +
               "is not associated with a STR. The maximum number of barcode " +
               "associated per STR is {maximum} and the minimum number " +
               "of barcode associated per STR is {minimum}. After " +
               "filtering STRs that do not have a minimum of {thres} " +
               "barcodes associated, {num_barcode} unique barcode is associated " +
               "with {num_STR} unique STR.\n")
    
    filt = copy.deepcopy(count_filtered_df)
    filt["STR"] = filt["barcode"].map(association)

    total_bc = len(filt)
    no_STR = filt.isnull().sum()["STR"]
    filt = filt.dropna()
    maximum = filt.STR.value_counts().max()
    minimum = filt.STR.value_counts().min()

    # filtering
    STR_count = filt.STR.value_counts().to_dict()
    filt["bc_count"] = filt["STR"].map(STR_count)
    filt = filt.fillna(0)
    group = filt[filt["bc_count"]>=min_barcode]

    # output matrix
    group.to_csv(out_dir + "rep" + str(group_num) + "_count_matrix.csv", index=False)

    group = group[group.columns.difference(['bc_count'])]

    num_bc = len(group)
    num_STR = len(group.STR.unique())

    # print the output messages 
    percent_not_associated = "{:.2f}".format((no_STR/total_bc)*100)
    
    print("Replicate " + str(group_num) + ":")
    print(out_msg.format(total=total_bc, no_STR=no_STR,
                         percent=percent_not_associated,
                         maximum=maximum, minimum=minimum,
                         thres=min_barcode, 
                         num_barcode=num_bc, num_STR=num_STR),
          flush=True)
    
    #record the statistic
    sum_dict["replicate" + str(group_num)]["associated_barcode"] = total_bc - no_STR
    sum_dict["replicate" + str(group_num)]["percent_associated_barcode"] = "{:.2f}".format(((total_bc - no_STR)/total_bc)*100)
    sum_dict["replicate" + str(group_num)]["fin_barcode"] = num_bc
    sum_dict["replicate" + str(group_num)]["fin_STR"] = num_STR
    
    return group
    
# plotting related function
def gDNA_scatter (unfilt_all, rep_num, out_dir):
    fig = plt.figure(figsize=((rep_num + 1)*10, 12))
    ax_num = 100 + rep_num * 10
    
    for i in range (1, rep_num):
        for j in range (i+1, rep_num + 1):
            ax_num += 1
            ax = sns.regplot(data=unfilt_all, x="gDNA"+str(i), y="gDNA"+str(j),
                             color="orange", ax=fig.add_subplot(ax_num));
            ax.set(xlabel='gDNA' + str(i) + ' barcode count',
                   ylabel='gDNA' + str(j) + ' barcode count')
            ax.set_title("Pearson R is " + str("{:.2f}".format(stats.pearsonr(unfilt_all["gDNA"+str(i)],
                                                               unfilt_all["gDNA"+str(j)])[0])));
    
    plt.savefig(out_dir + "gDNA_scatter.png")

def cDNA_scatter (unfilt_all, rep_num, out_dir):
    fig = plt.figure(figsize=((rep_num + 1)*10, 12))
    ax_num = 100 + rep_num * 10
    
    for i in range (1, rep_num):
        for j in range (i+1, rep_num + 1):
            ax_num += 1
            ax = sns.regplot(data=unfilt_all, x="cDNA"+str(i), y="cDNA"+str(j),
                             color="orange", ax=fig.add_subplot(ax_num));
            ax.set(xlabel='cDNA' + str(i) + ' barcode count',
                   ylabel='cDNA' + str(j) + ' barcode count')
            ax.set_title("Pearson R is " + str("{:.2f}".format(stats.pearsonr(unfilt_all["cDNA"+str(i)],
                                                               unfilt_all["cDNA"+str(j)])[0])));
    
    plt.savefig(out_dir + "cDNA_scatter.png")

def ratio_scatter_unfilt (dfs, rep_num, kind, out_dir):
    fig = plt.figure(figsize=((rep_num + 1)*10, 12))
    ax_num = 100 + rep_num*10
    shared = None

    for i in range(0, rep_num):
        df = dfs[i]
        cDNA = "cDNA" + str(i+1)
        gDNA = "gDNA" + str(i+1)

        ax_num += 1

        if not shared:
            ax = sns.regplot(data=df, x=gDNA, y=cDNA,
                             color="orange", 
                             ax=fig.add_subplot(ax_num));
            shared = ax
        else:
            ax = sns.regplot(data=df, x=gDNA, y=cDNA,
                             color="orange", 
                             ax=fig.add_subplot(ax_num, 
                                                sharex=shared,
                                                sharey=shared));
        ax.set(xlabel= gDNA + ' barcode count', ylabel= cDNA + ' barcode count')
        ax.set_title("Replicate " + str(i+1));
    
    plt.savefig(out_dir + kind + "_ratio.png")
    
def gDNA_distribution (unfilt_all, rep_num, out_dir):
    plt.figure()
    gDNA_cols = []
    for i in range(0, rep_num):
        gDNA_cols.append("gDNA" + str(i+1))

    ax = sns.histplot(data=unfilt_all[gDNA_cols],
                      multiple="dodge", shrink=0.8, discrete=True)
    ax.set_title("gDNA barcode count distribution - raw count");
    ax.set(xlabel= "gDNA barcode count");

    ax2 = plt.axes([0.24, 0.24, 0.55, 0.55])
    ax2 = sns.histplot(data=unfilt_all[gDNA_cols], ax=ax2,
                       multiple="dodge", shrink=0.8, discrete=True)
    ax2.set_title('zoom')
    ax2.set_xlim([0,150]);
    ax2.set_ylim([0,100000]);
    
    plt.savefig(out_dir + "gDNA_distribution.png")
    
def cDNA_distribution (unfilt_all, rep_num, out_dir):
    plt.figure()
    cDNA_cols = []
    for i in range(0, rep_num):
        cDNA_cols.append("cDNA" + str(i+1))

    ax = sns.histplot(data=unfilt_all[cDNA_cols],
                      multiple="dodge", shrink=0.8, discrete=True)
    ax.set_title("cDNA barcode count distribution - raw count");
    ax.set(xlabel= "cDNA barcode count");

    ax2 = plt.axes([0.24, 0.24, 0.55, 0.55])
    ax2 = sns.histplot(data=unfilt_all[cDNA_cols], ax=ax2,
                       multiple="dodge", shrink=0.8, discrete=True)
    ax2.set_title('zoom')
    ax2.set_xlim([0,500]);
    ax2.set_ylim([0,2000]);
    
    plt.savefig(out_dir + "cDNA_distribution.png")

def getargs():
    parser = argparse.ArgumentParser()
    
    # input file & output directory 
    inout_group = parser.add_argument_group("Input/Output")
    inout_group.add_argument("--seqdir", help="Directory of all the expression read sequencing result", 
                             type=str, required=True)
    inout_group.add_argument("--names", help="A separate .txt file that contains all the sample names", 
                             type=str, required=True)
    inout_group.add_argument("--filetype", help="File type of reads", type=str, 
                             choices=["fastq.gz", "fastq", "fq"],
                             required=True)
    inout_group.add_argument("--numreplicate", help="Number of replicates used", type=int,
                             required=True)
    inout_group.add_argument("--association", help="Path to association.tsv file", type=str,
                             required=True)
    inout_group.add_argument("--outdir", help="Path to output directory", type=str,
                             required=True)
    inout_group.add_argument("--export_fastq", help="Whether or not to output read file, default is No", 
                             type=str, choices=["Yes","No"], default="No")
    
    # levenshtein filter 
    filter_group = parser.add_argument_group("Filtering related threshold")
    filter_group.add_argument("--read_length",
                              help="Length of the expression read, default is 101 (NovaSeq)",
                              type=int, default=101)
    filter_group.add_argument("--lev_match",
                              help="Number of base pairs of exact match required after the 20bp barcode",
                              type=int, default=5)
    filter_group.add_argument("--lev_thres",
                              help="Maximum levenshtein score required to kept the read, default is 5",
                              type=int, default=5)
    filter_group.add_argument("--min_count",
                              help="Minimum read count required to kept the read, default is 10",
                              type=int, default=10)
    filter_group.add_argument("--min_barcode",
                              help="Minimum number of barcode required to be associated per STR, default is 3",
                              type=int, default=3)
    
    # get argument
    args = parser.parse_args()
    
    return args

def main(args):
    
    # input/output
    seq_dir = args.seqdir
    input_files = args.names
    file_type = args.filetype
    rep_num = args.numreplicate
    association_path = args.association
    out_dir = args.outdir
    if args.export_fastq == "Yes":
        export_file = True
    if args.export_fastq == "No":
        export_file = False

    # filtering related parameters 
    # read length 
    read_len = args.read_length
    # levenshtein
    lev_matching_length = args.lev_match
    lev_threshold = args.lev_thres
    # minimum read count
    min_count = args.min_count
    min_barcode = args.min_barcode
    
    # checking input file existence
    if not os.path.exists(input_files):
        common.WARNING("Error: %s does not exist"%input_files)
        return 1    

    # checking if association.tsv exists
    if not os.path.exists(association_path):
        common.WARNING("Error: %s does not exist"%association_path)
        return 1

    # checking if out_dir exists, if not, create the out_dir
    if not os.path.exists(os.path.dirname(out_dir)):
         os.mkdir(os.path.dirname(out_dir))
    
    # checking if plot_dir exists,
    # if not create /plot/ directory in the out_dir
    plot_dir = out_dir + "plot/"
    if not os.path.exists(os.path.dirname(plot_dir)):
        os.mkdir(os.path.dirname(plot_dir))

    # parsing the input file 
    cDNA_names = []
    gDNA_names = []
    file = open(input_files, "r")
    for line in file:

        if line.rstrip() == "#cDNA":
            reading_cDNA = True
            reading_gDNA = False
            continue
        if line.rstrip() == "#gDNA":
            reading_cDNA = False
            reading_gDNA = True
            continue

        # checking read file existences 
        if reading_cDNA:
            cDNA_path = seq_dir + line.rstrip()
            if not os.path.exists(cDNA_path):
                print("Error: %s does not exist"%cDNA_path)
                return 1
            else:
                cDNA_names.append(line.rstrip())
        if reading_gDNA:
            gDNA_path = seq_dir + line.rstrip()
            if not os.path.exists(gDNA_path):
                print("Error: %s does not exist"%gDNA_path)
                return 1
            else:
                gDNA_names.append(line.rstrip())

    # checking whether:
    # 1. cDNA list and gDNA list have matched length
    # 2. matches the input replicate number 
    if len(cDNA_names) != len(gDNA_names):
        print("Error: The amount of input cDNA files and gDNA files does not match")

        if len(cDNA_names) != rep_num:
            print("Error: The amount of input cDNA files does not match anticipated replicate number")

        if len(gDNA_names) != rep_num:
            print("Error: The amount of input gDNA files does not match anticipated replicate number")

        return 1

    # summary statistics dictionry 
    sum_dict = {}

    # pre-processing and count the barcode 
    count_dicts = []
    for i in range(0, rep_num):
        cDNA_name = cDNA_names[i]
        gDNA_name = gDNA_names[i]
        group_num = i + 1

        # filter cDNA
        c_path = seq_dir + cDNA_name
        c_bc_dict = filter_read (c_path, file_type,
                                 "cDNA", read_len,
                                 group_num, out_dir,
                                 lev_matching_length,
                                 lev_threshold,
                                 sum_dict, export_file)

        # filter gDNA
        g_path = seq_dir + gDNA_name
        g_bc_dict = filter_read (g_path, file_type,
                                 "gDNA", read_len,
                                 group_num, out_dir,
                                 lev_matching_length,
                                 lev_threshold,
                                 sum_dict, export_file)

        # check if the barcode dict is empty
        if not c_bc_dict:
            print(c_path + " generate an empty barcode dictionary", flush=True)

        if not g_bc_dict:
            print(g_path + " generate an empty barcode dictionary", flush=True)

        # generate count dictionary 
        # column name based on replicate number
        cDNA = "cDNA" + str(group_num)
        gDNA = "gDNA" + str(group_num)

        bc_count = {}
        for barcode in (set(c_bc_dict.keys()).union(set(g_bc_dict.keys()))):
            if barcode not in bc_count:
                bc_count[barcode] = {cDNA:0, gDNA:0}

            if barcode in c_bc_dict:
                bc_count[barcode][cDNA] = c_bc_dict[barcode]
            if barcode in g_bc_dict:
                bc_count[barcode][gDNA] = g_bc_dict[barcode]

        count_dicts.append(bc_count)

    # generate the unfiltered count matrix
    print("start generating the unfiltered count matrix...",
          flush=True)
    count_dfs = []
    for i in range(0, rep_num):
        count_dfs.append(pd.DataFrame.from_dict(count_dicts[i], orient="index").reset_index())
        count_dfs[i] = count_dfs[i].rename(columns={"index": "barcode"})

        count_dfs[i].to_csv(out_dir + "unfilt_count_matrix" + str(i+1) + ".csv",
                            index=False)  
        print("finish generating unfiltered count matrix " + str(i+1) + "...",
              flush=True)
    print("done\n", flush=True)

    # create the association dictionary 
    association_dict = {}

    file = open(association_path)
    for line in file:
        if line.rstrip():
            info = list(line.rstrip().split("\t"))
            barcode = info[0]
            STR = info[1]

            association_dict[barcode] = STR
        else:
            file.close()
            break 

    # filter the count matrix and associate BC with STR
    groups = []
    for i in range(0, rep_num):
        group_num = i + 1
        count_df = count_dfs[i]
        filt = count_filter(count_df, min_count)

        # record statistics 
        sum_dict["replicate"+str(group_num)] = {
            "init_cg_barcode": len(count_df),
            "count_filtered_barcode": len(filt)
        }

        groups.append(BC_STR_association (filt, out_dir, 
                                          association_dict, min_barcode,
                                          group_num, sum_dict))

    # barcode adds up - generate matrix used for active STR identification 
    agg_groups = []
    for i in range(0, rep_num):
        gDNA = "gDNA" + str(i+1)
        cDNA = "cDNA" + str(i+1)
        agg_group = groups[i].groupby(["STR"]).agg({gDNA:"sum", cDNA:"sum"})
        agg_group = agg_group.reset_index()
        agg_groups.append(agg_group)

    for i in range(0, rep_num):
        if i == 0:
            agg_all = pd.merge(left=agg_groups[0], right=agg_groups[1],
                                how="outer",
                                left_on="STR", right_on="STR")
        elif i == 1:
            continue
        elif i == rep_num:
            break
        else:
            agg_all = pd.merge(left=agg_all, right=agg_groups[i],
                                how="outer",
                                left_on="STR", right_on="STR")

    agg_all = agg_all.fillna(0)
    print("start generating the aggregated count matrix...")
    agg_all.to_csv(out_dir + "aggregate_count_matrix.csv", index=False)
    print("done\n", flush=True)

    print("A total of " + str(len(agg_all)) + " STR is captured " +
          "when combining the " + str(rep_num) + " replicates. \n",
          flush=True) 
    
    # output the summary statistics 
    sum_df = pd.DataFrame.from_dict(sum_dict, orient="index").reset_index()
    sum_df.to_csv(out_dir + "summary.csv", index=False)
    
    print("Finished generating the matricies. \n\n\n",
          flush=True)
    
    # plot
    print("Start plotting...\n", flush=True)
    for i in range(0, rep_num):
        if i == 0:
            unfilt_all = pd.merge(left=count_dfs[0], right=count_dfs[1],
                                  how="outer",
                                  left_on="barcode", right_on="barcode")
        elif i == 1:
            continue
        elif i == rep_num:
            break
        else:
            unfilt_all = pd.merge(left=unfilt_all, right=count_dfs[i],
                                  how="outer",
                                  left_on="barcode", right_on="barcode")

    unfilt_all = unfilt_all.fillna(0)
    
    gDNA_scatter (unfilt_all, rep_num, plot_dir)
    cDNA_scatter (unfilt_all, rep_num, plot_dir)
    ratio_scatter_unfilt (count_dfs, rep_num, "raw", plot_dir)
    gDNA_distribution(unfilt_all, rep_num, plot_dir)
    cDNA_distribution(unfilt_all, rep_num, plot_dir)
    ratio_scatter_unfilt (groups, rep_num, "associated", plot_dir)  
    
    print("Done.", flush=True)
    
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