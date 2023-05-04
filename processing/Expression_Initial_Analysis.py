#!/usr/bin/env python3
"""
script that perform initial analysis on the PROCESSED expression result, 
including identifying active STRs and correlation calculation.

currently support initial analysis on the hSTR library and Uber library. 
"""

# set up
import os
import sys 
import copy 
import gzip
import argparse
import subprocess
import numpy as np
import pandas as pd 
import scipy.stats as stats

# helper functions
def hSTR_STR_type(in_df):
    
    df = copy.deepcopy(in_df)
    STR_name = df["STR"].str.rsplit("_", 1, expand=True)
    
    return set(STR_name[0])

def hSTR_STR_split(in_df):
    
    df = copy.deepcopy(in_df)
    
    STR_name = df["STR"].str.rsplit("_", 1, expand=True)
    df["STR"] = STR_name[0]
    df["type"] = STR_name[1]
    
    return df

def is_active (STR, sig_STRs):
    if STR in sig_STRs:
        return "active"
    else:
        return "inactive"

def uber_STR_split (in_df):
    df = copy.deepcopy(in_df)

    STR_name = df["STR"].str.split("_", 3, expand=True)
    df["STR"] = (STR_name[0].astype(str) + "_" + 
                 STR_name[1].astype(str) + "_" + STR_name[2].astype(str))
    
    df["type"] = STR_name[3]
    
    return df

def uber_get_strand (repeat_type):
    repeat_info = list(repeat_type.split("_"))
    
    if len(repeat_info) == 2:
        return -1
    else:
        return int(repeat_info[3].strip("*"))

def get_repeat_frequency (repeat_type):
    if "p5" in repeat_type:
        return 5
    
    elif "m5" in repeat_type:
        return -5
    
    elif "p3" in repeat_type:
        return 3
    
    elif repeat_type.endswith("ref"):
        return 0
    
    else: 
        return int(repeat_type.split("_")[1])  

def uber_get_motif (repeat_type):
    
    repeat_info = list(repeat_type.split("_"))
    
    if len(repeat_info) == 2 or repeat_info[2] == "ref":
        return repeat_info[0]
    
    else: 
        return repeat_info[2]

def ratio_pearson_correlation (in_df, rep_int, mode, motif_dict=None, frequency_dict=None):    
    df = copy.deepcopy(in_df)
    
    col_cDNA = "cDNA" + str(rep_int)
    col_gDNA = "gDNA" + str(rep_int)

    calculate_dict = {}
    ratio_dict = {}
    correlation_dict = {}
    original_motif_dict = {}

    for barcode, cDNA, gDNA, STR, bc_count, type_STR in zip (df.barcode, 
                                                             df[col_cDNA],
                                                             df[col_gDNA],
                                                             df.STR, 
                                                             df.bc_count,
                                                             df.type):
        if mode == "STR":
            motif = motif_dict[STR + "_" + type_STR]
            original_motif = motif
            frequency = get_repeat_frequency(type_STR)
        if mode == "Uber":
            motif = uber_get_motif (type_STR)
            original_motif = type_STR.split("_")[0]
            
            if frequency_dict:
                if (STR + "_" + type_STR) in frequency_dict:
                    frequency = frequency_dict[STR + "_" + type_STR]
                else:
                    frequency = get_repeat_frequency(type_STR)
        
        original_motif_dict[STR] = original_motif
        
        col_ratio = "ratio" + str(rep_int)
        ratio = cDNA/gDNA
        
        if barcode not in ratio_dict:
            ratio_dict[barcode] = {
                "STR": "",
                "motif": "",
                "original motif": "",
                "repeat frequency": 0,
                col_ratio: 0
            }
        ratio_dict[barcode]["STR"] = STR
        ratio_dict[barcode]["motif"] = motif
        ratio_dict[barcode]["original motif"] = original_motif
        ratio_dict[barcode]["repeat frequency"] = frequency
        ratio_dict[barcode][col_ratio] = ratio
        
        key = STR + "_" + motif
        if key not in calculate_dict:
            calculate_dict[key] = {
                "ratios":[],
                "repeat_frequency":[]
            }
        calculate_dict[key]["ratios"].append(ratio)
        calculate_dict[key]["repeat_frequency"].append(frequency)

        if key not in correlation_dict:
            correlation_dict[key] = {
                "beta": 0,
                "pearson_r": 0,
                "pearson_pval": 0,
                "stderr": 0,
                "STR_type":set(),
                "bc_count": {},
            }

        correlation_dict[key]["STR_type"].add(type_STR)
        correlation_dict[key]["bc_count"][type_STR] = bc_count

    for key in calculate_dict:
        types = calculate_dict[key]["repeat_frequency"]
        ratios = calculate_dict[key]["ratios"]
        
        #use linear regression to caculate:
        #    1) beta
        #    2) pearson_r
        #    3) p-value
        #    4) standard error  
        stat_res = stats.linregress(x=types, y=ratios)
        correlation_dict[key]["beta"] = stat_res.slope
        correlation_dict[key]["pearson_r"] = stat_res.rvalue
        correlation_dict[key]["pearson_pval"] = stat_res.pvalue
        correlation_dict[key]["stderr"] = stat_res.stderr

    ratio_df = pd.DataFrame(ratio_dict).transpose()
    ratio_df = ratio_df.reset_index()
    ratio_df = ratio_df.rename(columns={'index':'barcode'})
    if mode == "hSTR":
        ratio_df = ratio_df[ratio_df.columns.drop(['original_motif'])]
    
    correlation_df = pd.DataFrame(correlation_dict).transpose()
    correlation_df = correlation_df.reset_index()
    correlation_df = correlation_df.rename(columns={'index':'STR'})
    STR_motif = correlation_df["STR"].str.rsplit("_", 1, expand=True)
    correlation_df["STR"] = STR_motif[0]
    correlation_df.insert(1, "motif", STR_motif[1])
    if mode == "Uber":
        correlation_df.insert(2, "original_motif",
                              correlation_df["STR"].map(original_motif_dict))
    
    return ratio_df, correlation_df 
def combine_statistics(betas, correlations, pvals, stderrs):
    
    num_ele = len(betas)
    for stat in (correlations, pvals, stderrs):
        if len(stat) != num_ele:
            print(stat + "have different numbers of elements")
            
            return 1
    
    # if only 1 value, return as it is 
    if num_ele == 1:
        return betas[0], stderrs[0], correlations[0], pvals[0], num_ele
            
    
    
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
#               "meta_p: " + str(final_p) + "\n")
    
    return final_beta, final_se, final_corr, final_p, num_ele

def combined_result (correlations, mode):
    statistics = {}
    for correlation in correlations:
        
        if mode == "hSTR":
            for STR, beta, corr, pval, stderr, motif in zip(correlation.STR, 
                                                            correlation.beta,
                                                            correlation.pearson_r,
                                                            correlation.pearson_pval,
                                                            correlation.stderr, 
                                                            correlation.motif):
                key = STR + "_" + motif

                if key not in statistics:
                    statistics[key] = {
                        "betas": [],
                        "corrs": [],
                        "pvals": [],
                        "stderrs": [],
                    }

                statistics[key]["betas"].append(beta)
                statistics[key]["corrs"].append(corr) 
                statistics[key]["pvals"].append(pval) 
                statistics[key]["stderrs"].append(stderr) 
                
        if mode == "Uber":
            for STR, beta, corr, pval, stderr, motif, ori, strand in zip(correlation.STR, 
                                                                         correlation.beta,
                                                                         correlation.pearson_r,
                                                                         correlation.pearson_pval,
                                                                         correlation.stderr, 
                                                                         correlation.motif,
                                                                         correlation.original_motif,
                                                                         correlation.strand):
                
                key = STR + "_" + motif + "_" + ori + "_" + str(strand)

                if key not in statistics:
                    statistics[key] = {
                        "betas": [],
                        "corrs": [],
                        "pvals": [],
                        "stderrs": [],
                    }

                statistics[key]["betas"].append(beta)
                statistics[key]["corrs"].append(corr) 
                statistics[key]["pvals"].append(pval) 
                statistics[key]["stderrs"].append(stderr) 
        
    combined_statistics = {
        "STR": [],
        "motif": [],
        "beta": [],
        "correlation": [],
        "p-value": [],
        "stderr": [],
        "num_rep": []
    }
    
    for key in statistics:
        STR_motif = key.split("_")
        STR = STR_motif[0] + STR_motif[1] + STR_motif[2]
        motif = STR_motif[3]

        betas = statistics[key]["betas"]
        corrs = statistics[key]["corrs"]
        pvals = statistics[key]["pvals"]
        stderrs = statistics[key]["stderrs"]

        # combined the statistics
        final_beta, final_se, final_corr, final_p, num_ele = combine_statistics(betas, corrs, pvals, stderrs)

        combined_statistics["STR"].append(STR)
        combined_statistics["motif"].append(motif)
        combined_statistics["beta"].append(final_beta)
        combined_statistics["correlation"].append(final_corr)
        combined_statistics["p-value"].append(final_p)
        combined_statistics["stderr"].append(final_se)
        combined_statistics["num_rep"].append(num_ele)        
        
        if mode == "Uber":
            if "original motif" not in combined_statistics:
                combined_statistics["original motif"] = []
            
            if "strand" not in combined_statistics:
                combined_statistics["strand"] = []
            
            original_motif = STR_motif[4]
            combined_statistics["original motif"].append(original_motif)
            
            strand = STR_motif[5]
            combined_statistics["strand"].append(strand)
                
    stat_df = pd.DataFrame(combined_statistics)

    return stat_df

def getargs():
    parser = argparse.ArgumentParser()
    
    # processing mode 
    mode_group = parser.add_argument_group("Mode")
    mode_group.add_argument("--mode", help="Processing mode",
                            choices=["hSTR", "Uber"],
                            type=str, required=True)
    
    # input file & output directory 
    inout_group = parser.add_argument_group("Input/Output")
    inout_group.add_argument("--datadir", help="Directory of all the expression count matricies", 
                             type=str, required=True)
    inout_group.add_argument("--numreplicate", help="Number of replicates used", type=int,
                             required=True)
    inout_group.add_argument("--DESeq2", help="Path to MPRA_DESeq2.r file, required if mode=hSTR",
                             type=str, required=False)
    inout_group.add_argument("--refSTR", help="Path to Mikhail's array_probes_split.tsv file, required if mode=hSTR",
                             type=str, required=False)
    inout_group.add_argument("--refFrequency", help="Path to actual repeat frequency for probe selected from hSTR library, required if mode=Uber",
                             type=str, required=False)
    inout_group.add_argument("--outdir", help="Path to output directory", type=str,
                             required=True)
    
    # get argument
    args = parser.parse_args()
    
    return args

def main(args):
    # parameters 
    # processing mode 
    mode = args.mode 
    
    # input/output 
    data_dir = args.datadir
    rep_num = args.numreplicate
    out_dir = args.outdir
    
    # checking if data_dir exists
    if not os.path.exists(os.path.dirname(data_dir)):
        print("Error: data directory %s does not exists"%data_dir)
        return 1
    
    # checking if out_dir exists, if not, create the out_dir
    if not os.path.exists(os.path.dirname(out_dir)):
        os.mkdir(os.path.dirname(out_dir))
    
    # checking if the mode-specific required file exists:
    if mode == "hSTR":
        DESeq2_script_path = args.DESeq2
        ref_info_path = args.refSTR
        
        if not os.path.exists(DESeq2_script_path):
            print("Error: %s does not exist"%DESeq2_script_path)
            return 1
        
        if not os.path.exists(ref_info_path):
            print("Error: %s does not exist"%ref_info_path)
            return 1
        
        # DESeq2
        print("mode is hSTR...\n", flush = True)
        print("start performing DESeq2...", flush=True)
        DESeq2_args = [data_dir + "aggregate_count_matrix.csv",
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

        sig_STR_df = deseq2[deseq2["pvalue"] <= 0.01/len(hSTR_STR_type(deseq2))]
        sig_full_STRs = list(sig_STR_df.STR)

        sig_STR_df = hSTR_STR_split(sig_STR_df)
        sig_STRs = set(sig_STR_df.STR)

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

        characterization = pd.read_csv(normalized_count)
        characterization = characterization.rename(columns={"Unnamed: 0": "STR"})
        characterization["active"] = characterization["STR"].apply(is_active, args=([sig_full_STRs]))
        characterization["motif"] = characterization.STR.map(STR_motif)
        characterization["str_len"] = characterization.STR.map(STR_len)
        characterization["str_len"] = characterization["str_len"].fillna(0)
        characterization["str_len"] = characterization["str_len"].astype(int)

        print("start generating characterization.csv file based on reporter activity...",
              flush=True)
        characterization.to_csv(out_dir + "characterization" + ".csv",
                                index=None) 
        print("done\n", flush=True)        
    
    if mode == "Uber":
        ref_frequency_path = args.refFrequency
        
        if not os.path.exists(ref_frequency_path):
            print("Error: %s does not exist"%ref_frequency_path)
            return 1
        
        #probe frequency reference 
        freq_ref = pd.read_csv(ref_frequency_path)
        freq_dict = dict((zip(freq_ref.STR, freq_ref.repeat_freq)))
    
    # pearson correlation calculation and matrix generation 
    print("start calculating perason correlation of expression and STR type...",
          flush=True)
    ratio_dfs = []
    correlation_dfs = []
    for i in range(0, rep_num):
        bc_group = pd.read_csv(data_dir + "rep" + str(i+1) + "_count_matrix.csv")
    
        if mode == "hSTR":
            bc_group = hSTR_STR_split(bc_group)

            # select significant STRs
            bc_group = bc_group[bc_group.STR.isin(sig_STRs)]            
            
            ratio_df, correlation_df = ratio_pearson_correlation(bc_group, i+1, mode,
                                                                 motif_dict=STR_motif)
            ratio_dfs.append(ratio_df)
            correlation_dfs.append(correlation_df)

        if mode == "Uber":
            bc_group = uber_STR_split(bc_group)
            bc_group["strand"] = bc_group["type"].apply(uber_get_strand)

            ratio_ori, correlation_ori = ratio_pearson_correlation(bc_group[bc_group["strand"] == -1],
                                                                   i+1, mode, frequency_dict=freq_dict)
            ratio_ori["strand"] = -1
            correlation_ori["strand"] = -1

            ratio_0, correlation_0 = ratio_pearson_correlation(bc_group[bc_group["strand"] == 0],
                                                               i+1, mode, frequency_dict=freq_dict)
            ratio_0["strand"] = 0
            correlation_0["strand"] = 0

            ratio_1, correlation_1 = ratio_pearson_correlation(bc_group[bc_group["strand"] == 1],
                                                               i+1, mode, frequency_dict=freq_dict)
            ratio_1["strand"] = 1
            correlation_1["strand"] = 1

            ratio_df = pd.concat([ratio_ori, ratio_0, ratio_1])
            correlation_df = pd.concat([correlation_ori, correlation_0, correlation_1])

            ratio_dfs.append(ratio_df)
            correlation_dfs.append(correlation_df)
       
        # output matrix
        print("start generating final ratio matrix and correlation matrix " + str(i+1) + "...",
              flush=True)
        ratio_df.to_csv(out_dir + "rep" + str(i+1) + "_ratio_matrix.csv", index=None)
        correlation_df.to_csv(out_dir + "rep" + str(i+1) + "_pearson_correlation_matrix.csv",
                              index=None)
        print("done\n", flush=True)
        
    # combined the replicates
    print("start combining the matricies ...",
          flush=True)
    combined_df = combined_result(correlation_dfs, mode)
    combined_df.to_csv(out_dir + "combined_results.csv", index=None)
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







    
    
    