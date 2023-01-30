#!/usr/bin/env python3
"""
script for processing the expression reads
which will output multiple types of count matrix 
to the desinated output directory 
"""

# set up
import os
import sys
import copy
import gzip
import argparse
import subprocess
import Levenshtein
import numpy as np
import pandas as pd
import scipy.stats as stats

# helper functions
def process(lines=None):
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
    
    # the full identical sequence that could be 
    # used for qc
    qc_seq = "TCTAGAGGTTCGTCGACGCGATTATTATCATTACTTGTACAGCTCGTCCATGCCGAGAGTGATC"
    
    # the sqeuence required to be exact match 
    exact_match = qc_seq[:levenshtein_matching_length]
    
    # the sequence being checked 
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
                 lev_threshold):
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
            # output file
            file_out = gzip.open(out_dir + fname, "wb")
            file_fail = gzip.open(out_dir + fail, "wb")
        else:
            # input file
            file_in = open(path, "r")
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
                    
                    # output read that fail the qc
                    if ".gz" in file_type:
                        file_fail.write((read_id + "\n" + seq + "\n" + 
                                         sign + "\n" + qual + "\n").encode())
                    else:
                        file_fail.write(read_id + "\n" + seq + "\n" + 
                                        sign + "\n" + qual + "\n")
                
                # empty the line list 
                lines = []
        
        # close the files
        file_in.close()
        file_out.close()
        file_fail.close()
        
        # print output message
        print(txt.format(filt = filtered_read,
                         total = total_read, 
                         R = read_type,
                         lev_len = lev_matching_length,
                         lev_thres = lev_threshold,
                         remain = remained_read,
                         percent = "{:.2f}".format((remained_read/total_read)*100)),
              flush=True)
        print("finished writing filtered reads to " + out_dir + fname + "\n",
              flush=True)
        
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
                        group_num):

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

    print("Replicate " + str(group_num) + ":")
    print(out_msg.format(total=total_bc, no_STR=no_STR,
                         percent=str("{:.2f}".format((no_STR/total_bc)*100)),
                         maximum=maximum, minimum=minimum,
                         thres=min_barcode, 
                         num_barcode=num_bc, num_STR=num_STR),
          flush=True)
    
    return group

def STR_type(in_df):
    
    df = copy.deepcopy(in_df)
    STR_name = df["STR"].str.rsplit("_", 1, expand=True)
    
    return len(STR_name[0].unique())
    

def STR_split(in_df):
    
    df = copy.deepcopy(in_df)
    
    STR_name = df["STR"].str.rsplit("_", 1, expand=True)
    df["STR"] = STR_name[0]
    df["type"] = STR_name[1]
    
    return df

def calculate_ratio (df, rep_num):
    cDNA = "cDNA" + str(rep_num)
    gDNA = "gDNA" + str(rep_num)
    ratio = "ratio" + str(rep_num)
    
    df[ratio] = df[cDNA]/df[gDNA]

def is_active (STR, sig_STRs):
    if STR in sig_STRs:
        return "active"
    else:
        return "inactive"

def motif_sep (seq_motif):
    only_CG = True
    only_AT = True
    if "A" in seq_motif or "T" in seq_motif:
        only_CG = False        
    
    if "C" in seq_motif or "G" in seq_motif:
        only_AT = False
    
    if (not only_CG) and (not only_AT):
        return "mixed"
    else:
#         if only_CG:
#             return "CG_only"
#         if only_AT:
#             return "AT_only"
        return "CG or AT"

def type_to_int (type_str):
    if type_str == "ref":
        return int(0)
    else:
        if type_str[0] == "p":
            return int(type_str[1:])
        if type_str[0] == "m":
            return -int(type_str[1:])

def ratio_pearson_correlation (in_df, rep_int, motif_dict, str_len_dict):
    
    col_cDNA = "cDNA" + str(rep_int)
    col_gDNA = "gDNA" + str(rep_int)

    calculate_dict = {}
    out_dict = {}

    for barcode, cDNA, gDNA, STR, bc_count, type_STR in zip(in_df.barcode,
                                                            in_df[col_cDNA],
                                                            in_df[col_gDNA],
                                                            in_df.STR, 
                                                            in_df.bc_count,
                                                            in_df.type):
        if STR not in calculate_dict:
            calculate_dict[STR] = {
                "ratios":[],
                "type_int":[]
            }
        calculate_dict[STR]["ratios"].append(cDNA/gDNA)
        calculate_dict[STR]["type_int"].append(type_to_int(type_STR))

        if STR not in out_dict:
            out_dict[STR] = {
                "beta": 0,
                "pearson_r": 0,
                "pearson_pval": 0,
                "stderr": 0,
                "motif": "",
                "STR_type":set(),
                "bc_count": {},
                "str_len": {},
            }
        out_dict[STR]["STR_type"].add(type_STR)
        out_dict[STR]["bc_count"][type_STR] = bc_count
        out_dict[STR]["motif"] = motif_dict[STR + "_" + type_STR]
        out_dict[STR]["str_len"][type_STR] = str_len_dict[STR + "_" + type_STR]


    for STR in calculate_dict:
        types = calculate_dict[STR]["type_int"]
        ratios = calculate_dict[STR]["ratios"]
         
        #directly calculated pearson_r
        #out_dict[STR]["pearson_r"] = stats.pearsonr(x=types, y=ratios)[0]
        #out_dict[STR]["pearson_pval"] = stats.pearsonr(x=types, y=ratios)[1]
        
        #use linear regression to caculate:
        #    1) beta
        #    2) pearson_r
        #    3) p-value
        #    4) standard error  
        stat_res = stats.linregress(x=types, y=ratios)
        out_dict[STR]["beta"] = stat_res.slope
        out_dict[STR]["pearson_r"] = stat_res.rvalue
        out_dict[STR]["pearson_pval"] = stat_res.pvalue
        out_dict[STR]["stderr"] = stat_res.stderr

    # select for ones that are in groups of 3 or more
    for STR in list(out_dict.keys()):
        if len(out_dict[STR]["STR_type"]) < 3:
            out_dict.pop(STR)

    #display(out_dict)
            
    out_df = pd.DataFrame(out_dict).transpose()
    out_df = out_df.reset_index()
    out_df = out_df.rename(columns={'index':'STR'})
    
    return out_df

def combine_statistics(betas, correlations, pvals, stderrs):
    
    num_ele = len(betas)
    for stat in (correlations, pvals, stderrs):
        if len(stat) != num_ele:
            print(stat + "have different numbers of elements")
    
    ws = []
    bws = []
    cws = []
    
    for i in range(0, num_ele):
        w = 1/(stderrs[i] ** 2)
        
        ws.append(w)
        bws.append(betas[i] * w)
        cws.append(correlations[i] * w)
      
    final_se = (1/sum(ws)) ** 0.5
    final_beta = sum(bws)/sum(ws)
    final_corr = sum(cws)/sum(ws)
    
    # use the above meta-analysis formula
    z = final_beta/final_se
    final_p = stats.norm.sf(abs(z))*2
    
#     # check
#     if num_ele == 1:
#         print("beta: " + str(final_beta) + "\n" +
#               "correlation: " + str(final_corr) + "\n" + 
#               "stderr: " + str(final_se) + "\n" +
#               "meta_p: " + str(final_p) + "\n" + 
#               "scipy_p: " + str(combined_p) + "\n")
    
    return final_beta, final_se, final_corr, final_p, num_ele

def stat_data (expressions, rep_num):
    statistics = {}
    for i in range(0, rep_num):
        expression = expressions[i]


        for STR, beta, corr, pval, stderr, motif, ID, chrom, pos, end in zip(expression.STR, 
                                                                             expression.beta, 
                                                                             expression.pearson_r, 
                                                                             expression.pearson_pval, 
                                                                             expression.stderr, 
                                                                             expression.motif,
                                                                             expression.gene_id, 
                                                                             expression.str_chr, 
                                                                             expression.str_pos,
                                                                             expression.str_end):

            if STR not in statistics:
                statistics[STR] = {
                    "betas": [],
                    "corrs": [],
                    "pvals": [],
                    "stderrs": [],
                    "motif": "",
                    "gene_id": "",
                    "str_chr": "",
                    "str_pos": "",
                    "str_end": ""
                }

            statistics[STR]["betas"].append(beta)
            statistics[STR]["corrs"].append(corr) 
            statistics[STR]["pvals"].append(pval) 
            statistics[STR]["stderrs"].append(stderr) 
            statistics[STR]["motif"] = motif
            statistics[STR]["gene_id"] = ID
            statistics[STR]["str_chr"] = chrom
            statistics[STR]["str_pos"] = pos
            statistics[STR]["str_end"] = end

    combined_statistics = {
        "STR": [],
        "beta": [],
        "correlation": [],
        "p-value": [],
        "stderr": [],
        "num_rep": [],
        "motif": [],
        "gene_id": [],
        "chrom": [],
        "str.start": [],
        "str.end": []
    }

    for STR in statistics:
        betas = statistics[STR]["betas"]
        corrs = statistics[STR]["corrs"]
        pvals = statistics[STR]["pvals"]
        stderrs = statistics[STR]["stderrs"]

        # combined the statistics
        final_beta, final_se, final_corr, final_p, num_ele = combine_statistics(betas, corrs, pvals, stderrs)

        combined_statistics["STR"].append(STR)
        combined_statistics["beta"].append(final_beta)
        combined_statistics["correlation"].append(final_corr)
        combined_statistics["p-value"].append(final_p)
        combined_statistics["stderr"].append(final_se)
        combined_statistics["num_rep"].append(num_ele)

        combined_statistics["motif"].append(statistics[STR]["motif"])
        combined_statistics["gene_id"].append(statistics[STR]["gene_id"])
        combined_statistics["chrom"].append(statistics[STR]["str_chr"])
        combined_statistics["str.start"].append(statistics[STR]["str_pos"])
        combined_statistics["str.end"].append(statistics[STR]["str_end"])

    stat_df = pd.DataFrame(combined_statistics)

    return stat_df

def getargs():
    parser = argparse.ArgumentParser()
    
    # processing mode 
    mode_group = parser.add_argument_group("Mode")
    mode_group.add_argument("--mode", help="Processing mode",
                            choices=["ALL", "QC"],
                            type=str, required=True)
    
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
    inout_group.add_argument("--DESeq2", help="Path to MPRA_DESeq2.r file, required if mode=ALL",
                             type=str, required=False)
    inout_group.add_argument("--refSTR", help="Path to Mikhail's array_probes_split.tsv file, required if mode=ALL",
                             type=str, required=False)
    inout_group.add_argument("--tss_str", help="Path to Mikhail's tss_str_pairs.tsv file, required if mode=ALL",
                             type=str, required=False)
    inout_group.add_argument("--ens_gene", help="Path to gene annotation ens_genes.tab file, required if mode=ALL",
                             type=str, required=False)
    inout_group.add_argument("--outdir", help="Path to output directory", type=str,
                             required=True)
    
    # levenshtein filter 
    filter_group = parser.add_argument_group("Filtering related threshold")
    filter_group.add_argument("--read_length",
                              help="Length of the expression read, default is 84",
                              type=int, default=84)
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
    # parameters 
    # processing mode 
    mode = args.mode
    
    # input/output
    seq_dir = args.seqdir
    input_files = args.names
    file_type = args.filetype
    rep_num = args.numreplicate
    association_path = args.association
    out_dir = args.outdir

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
    
    # if processing mode is ALL
    if mode == "ALL":
        DESeq2_script_path = args.DESeq2
        ref_info_path = args.refSTR
        tss_str_path = args.tss_str
        ens_gene_path = args.ens_gene
        other_paths = [
            DESeq2_script_path,
            ref_info_path, 
            tss_str_path, 
            ens_gene_path
        ]
        
        # checking if the above files exists
        for other_file in other_paths:
            if not os.path.exists(other_file):
                common.WARNING("Error: %s does not exist"%other_file)
                return 1

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
                                 lev_threshold)

        # filter gDNA
        g_path = seq_dir + gDNA_name
        g_bc_dict = filter_read (g_path, file_type,
                                 "gDNA", read_len,
                                 group_num, out_dir,
                                 lev_matching_length,
                                 lev_threshold)

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
        filt = count_filter(count_dfs[i], min_count)
        groups.append(BC_STR_association (filt, out_dir, 
                                          association_dict, min_barcode,
                                          group_num))
    # only continue the following process if mode=ALL
    
    if mode=="ALL":
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

        # DESeq2
        print("start performing DESeq2...", flush=True)
        DESeq2_args = [out_dir + "aggregate_count_matrix.csv",
                       out_dir, 
                       str(rep_num)]
        cmd = [DESeq2_script_path] + DESeq2_args
        subprocess.call(cmd)
        print("done\n", flush=True)

        # active STR identification
        deseq2_path = out_dir + "deseq2_result.csv"
        normalized_count = out_dir + "normalized_aggregate_count_matrix.csv"

        deseq2 = pd.read_csv(deseq2_path)
        deseq2 = deseq2.rename(columns={"Unnamed: 0": "STR"})
        norm_mc = pd.read_csv(normalized_count)
        norm_mc = norm_mc.rename(columns={"Unnamed: 0": "STR"})

        sig_STR_df = deseq2[deseq2["pvalue"] <= 0.01/STR_type(norm_mc)]
        sig_full_STRs = list(sig_STR_df.STR)

        sig_STR_df = STR_split(sig_STR_df)
        sig_STRs = list(sig_STR_df.STR)

        # generate active STR characterization realated .csv file
        ref_STRs = pd.read_csv(ref_info_path, sep="\t")
        human_STR = ref_STRs[ref_STRs["organism"]=="hg38"]
        human_STR["STR"] = human_STR["id"] + "_" + human_STR["allele"]
        human_STR = human_STR[["id", "allele", "STR", "motif", "str_seq"]]
        human_STR.columns = ["STR", "type", "full_STR", "motif", "str_seq"]
        human_STR["str_len"] = human_STR.str_seq.str.len()
        human_STR["str_len"] = human_STR["str_len"].fillna(0)

        STR_motif = dict(zip(human_STR.full_STR, human_STR.motif))
        STR_len = dict(zip(human_STR.full_STR, human_STR.str_len.astype(int)))

        characterization = copy.deepcopy(norm_mc)
        characterization["active"] = characterization["STR"].apply(is_active, args=([sig_full_STRs]))
        characterization["motif"] = characterization.STR.map(STR_motif)
        characterization["motif_group"] = characterization["motif"].apply(motif_sep)
        characterization["str_len"] = characterization.STR.map(STR_len)
        characterization["str_len"] = characterization["str_len"].fillna(0)
        characterization["str_len"] = characterization["str_len"].astype(int)

        print("start generating characterization.csv file based on reporter activity...",
              flush=True)
        characterization.to_csv(out_dir + "characterization" + ".csv",
                                index=False) 
        print("done\n", flush=True)

        # obtain TSS and other gene info
        tss_str_info = pd.read_csv(tss_str_path, sep="\t")
        ens_gene_dict = {}
        ens_gene_file = open(ens_gene_path, "r")
        for line in ens_gene_file:
            if line.rstrip():
                line_list = list(line.strip().split("\t"))
                gene_id = line_list[0]
                gene_name = line_list[1]
                ens_gene_dict[gene_id] = gene_name
            else:
                ens_gene_file.close()

        # pearson correlation calculation and matrix generation
        print("start calculating perason correlation of expression and STR type...",
              flush=True)
        expressions = []
        for i in range(0, rep_num):
            bc_group = pd.read_csv(out_dir + "rep" + str(i+1) + "_count_matrix.csv")
            bc_group = STR_split(bc_group)

            # select significant STRs
            bc_group = bc_group[bc_group.STR.isin(sig_STRs)]

            # pearson correlation calculation
            expression = ratio_pearson_correlation(bc_group, i+1, STR_motif, STR_len)
            expression = expression.merge(tss_str_info,
                                          left_on="STR", right_on="str_id", how = "left")
            expression = expression[expression.columns.drop(['str_id', 'organism'])]
            expression.insert(1, "gene", 
                              expression["gene_id"].map(ens_gene_dict))
            expressions.append(expression)
            
            # output matrix
            print("start generating final correlation matrix " + str(i+1) + "...",
                  flush=True)
            expression.to_csv(out_dir + "rep" + str(i+1) + "_pearson_correlation_matrix.csv")
        print("done\n", flush=True)
        
        # combined the replicates
        print("start combining the matricies ...",
              flush=True)
        effect_df = stat_data(expressions, rep_num)
        effect_df.to_csv(out_dir + "combined_results.csv", index="False")
        print("done\n", flush=True)
        
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











