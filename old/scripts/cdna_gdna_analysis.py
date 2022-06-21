import os
import re
import gzip
import argparse
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
from analyze_bams import rev_comp


TAGS="/storage/mlamkin/projects/eSTR-data/first_round_tags.txt"
FIG_PATH='/storage/mlamkin/projects/eSTR-MPRA-analysis/figures/cdna_gdna_analysis/'

UMI_MARKER="[ACGT]{7}AGTCGTCGTGCTGGAACA" # the sequence after the umi in the beginning to indicate its a proper read
TAG_MARKER="[ACGT]{10}TCTAGAATTATT"       # the sequence after the tag in the beginning of reverse read

# UMI READ: [N]{7} + AGTCGTCGTGCTGGAACA
# TAG READ: [N]{10} + TCTAGAATTATT
# TAG = [N]{10}
# UMI = [N]{7}
# Do a check on the specified tag and umi read to see if it matches with the pattern specified
# then collect the information from forward and reverse reads into separate files
# reverse contains tag and is crucial for collecting information forward is not necessary
def collect_reads(reads, markers, lengths, outputs):
    """
    params
    Note all of the below parameters should be in the format (tag, umi) when it comes to respective files
    reads - list which contains single (one file) or paired (two files, tag file first) which are gzipped
    markers - list of markers to filter out bad reads from each respective fastq file in reads
    lengths - tuple of lengths (tag length, umi length) to isolate from each read
    outputs - list of output file paths to write filtered reads
    """
    assert len(reads) == len(outputs)

    # open fastq files
    for ind, read_file in enumerate(reads):
        reads[ind] = gzip.open(read_file, 'r')

    # remove previous files with name
    for output in outputs:
        if os.path.isfile(output):
            os.remove(output)

    # open output files
    for ind, read_file in enumerate(outputs):
        outputs[ind] = open(read_file, 'a')

    tag_loc = 0           # store where the tag was found to be located so we can also filter by Phred Score at the same location
    filtered_reads = 0    # store the amount of reads that are present after the filtering process
    count = 0             # count to tell when to grab UMI and tag (May need to use regex for tag)
    skip_read = False     # determine whether its a bad read or good

    for bin_reads in zip(*reads):
        index = count % 4

        # list to store new reads (will be length 4)
        if index == 0:
            locs = [-1]*len(reads)     # [tag start] or [tag start, umi start]
            skip_read = False
            new_reads = [['']*4 for i in range(len(reads))]

        # found a bad read that doesnt match expected
        if skip_read:
            count += 1
            continue

        # decode reads because read as binary and store in new reads
        for ind, bin_read in enumerate(bin_reads):
            new_reads[ind][index] = bin_read.decode('ascii')

        # if the current line is the sequence grab umi from forwards read, tag from reverse read, and add to header
        # grab reverse read to check for SNPs and errors
        if index == 1:
            passed_marker = True
            seqs = [new[index] for new in new_reads]
            for ind, seq, marker in zip(range(len(reads)), seqs, markers):
                matched = re.search(marker, seq)
                if matched: locs[ind] = matched.start() 
                else: passed_marker = False

            if passed_marker:
                tag_umi = []
                # collect tag and umi
                for ind, seq, marker, feat_len in zip(range(len(reads)), seqs, markers, lengths):
                    tag_umi.append(seq[locs[ind]:locs[ind]+feat_len])

                # insert tag and umi into header
                for ind, header in enumerate([new[0] for new in new_reads]):
                    upd_header = header.rstrip('\n').replace(' ', '_') + '_' + '_'.join(tag_umi) + '\n'
                    new_reads[ind][0] = upd_header
            else:
                skip_read = True

        # filter tag by phred score with mim 30
        elif index == 3:
            good_phred = True
            for char in new_reads[0][index][locs[0]:locs[0]+lengths[0]]:
                if ord(char) < 63:
                    good_phred = False
                    break

            # output new reads to paired end fastq file if filter is passed
            if good_phred:
                filtered_reads += 1
                for ind, output in enumerate(outputs):
                    for row in new_reads[ind]:
                        output.write(row)
        count += 1

    total_reads = count/4
    usable_rate = filtered_reads * 1.0/total_reads
    print("Initial Read Count: %i reads"%(total_reads))
    print("Post Phred and Forward/Reverse Marker Filters Read Count: %i reads"%(filtered_reads))
    print("Rate of usable reads: %f"%(usable_rate))
    return


# go through the header of each read and grab the umi and header to generate stats from
def grab_umi_tag(read_file, tags, umi=False):
    count = 0
    tag_umi_dict = defaultdict(list)
    for line in read_file:
        index = count % 4
        if index == 0:
            if umi:
                tag_umi = line.rstrip('\n').split('_')[-2:]
                # need to take reverse complement of tag since on reverse strand
                tag_umi[0] = rev_comp(tag_umi[0])
                if tags.get(umi_tag[0], False):
                    tag_umi_dict['UMI'].append(tag_umi[1])
                    tag_umi_dict['Tag'].append(tag_umi[0])
            else:
                tag = rev_comp(line.rstrip('\n').split('_')[-1])
                if tags.get(tag, False):
                    tag_umi_dict['Tag'].append(tag)
        count += 1
    return pd.DataFrame(tag_umi_dict)


# Sort and collapse 
def collapse_umis(dna):
    pre_collapse_reads = dna.shape[0]
    dna.sort_values(["UMI", "Tag"], inplace=True)
    dna.drop_duplicates(["UMI", "Tag"], inplace=True)
    dna.sort_values(["UMI", "Tag"], inplace=True)
    post_collapse_reads = dna.shape[0]
    print("Number reads collapse on UMI:", pre_collapse_reads - post_collapse_reads)
    return dna


def normalize_tag_comparison(cdna_gdna, fig_path, figname, figures=False):
    def graph(counts_mtx, cols, title, output):
        fig, ax = plt.subplots(figsize=(7,7))

        # Filter out all zeros
        nonzero_loc = ~(counts_mtx == 0)
        nonzero_first = nonzero_loc[:,cols[0]]
        nonzero_second = nonzero_loc[:,cols[1]]
        nonzero_locs = np.logical_and(nonzero_first, nonzero_second)
        
        # get first and second tags to compare in graph
        first_tag = counts_mtx[nonzero_locs, cols[0]]
        second_tag = counts_mtx[nonzero_locs, cols[1]]
        slope, intercept, r_value, p_value, std_err = stats.linregress(first_tag, second_tag)

        # graph the scatter plot
        ax.scatter(first_tag, second_tag)
        ax.plot(first_tag, slope*first_tag + intercept, color='red')
        plt.text(1, -1, 'R^2 = %f'%(r_value**2))
        ax.set_xlabel("Tag %i"%(cols[0]+1))
        ax.set_ylabel("Tag %i"%(cols[1]+1))
        ax.set_title(title)
        #ax.set_xlim(left=0, right=4)
        #ax.set_ylim(bottom=0, top=4)
        plt.show()
        fig.tight_layout()
        fig.savefig(output, bbox_inches='tight')
        return

    # create numpy matrix of counts and take log of all values
    normalized_counts = _count_mtx(cdna_gdna, 'cDNA/gDNA')
    normalized_counts[normalized_counts == 0] += 1
    normalized_counts = np.log(normalized_counts)
   
    # create graphs that has counts for each tag of a particular allele
    if figures:
        for pair in [(0, 1), (0, 2), (1, 2)]:
            output_path = fig_path + 'cDNA_gDNA_count%ivs%i_correlation_%s.png'%(pair[0]+1, pair[1]+1, figname)
            graph(normalized_counts, pair, "Normalized Transcription Rate Correlation", output_path)
    return


# calculate variances within the 3 tag replicates
# remove everything with 0s
# if random randomly shuffle all elements then calc variance
def variance(cdna_gdna, rand=False):
    counts_mtx = _count_mtx(cdna_gdna, 'cDNA/gDNA')

    # remove all rows that contain at least one zero
    counts_mtx = counts_mtx[~np.any(counts_mtx == 0, axis=1),:]

    # additional filter to remove outliers
    counts_mtx = counts_mtx[~np.any(counts_mtx > 4, axis=1),:]

    # output the number of tag trios left post filter
    print("Number of tag replicates with no zeros:", counts_mtx.shape[0])

    # randomly shuffle the matrix to compare against original variances
    if rand:
        origin_shape = counts_mtx.shape
        flat_mtx = counts_mtx.flatten()
        np.random.shuffle(flat_mtx)
        counts_mtx = flat_mtx.reshape(origin_shape[0], origin_shape[1])

    # create array of variances based on columns of replicates
    counts_var = np.var(counts_mtx, axis=1)
    return counts_var


# create matrix of counts from cDNA/gDNA data
def _count_mtx(cdna_gdna, counts_col):
    # generate indices for each tag in a tuple (row, col) read in from file
    allele_file = open("/storage/mlamkin/projects/eSTR-data/finalized_first_round_mpra_STR_v2.0_updated_labels.txt", 'r')
    allele_file.readline()
    alleles = [[line.rstrip('\n').split(' ')[5],
                line.rstrip('\n').split('\t')[6],
                line.rstrip('\n').split('\t')[7]] for line in allele_file]
    
    tag_to_ind = {}
    ind_to_tag = {}
    
    count = 0
    for allele in alleles:
        row = count//3
        col = count%3
        tag_to_ind[allele[0]] = (row, col)
        ind_to_tag[row] = allele
        count += 1
   
    # make a 3xN matrix that stores normalized cDNA/gDNA counts for each individual tag
    counts_mtx = np.zeros((int(len(tag_to_ind.keys())/3), 3)) 
    for index, df_row in cdna_gdna.iterrows():
        row, col = tag_to_ind.get(df_row['Tag'], (None, None))
        if row == None and col == None:
            continue
        counts_mtx[row, col] = df_row[counts_col]
    return counts_mtx


def main():
    parser = argparse.ArgumentParser(description="Analyze read counts of specific gDNA and cDNA tags within fastq files.")
    parser.add_argument('-p', '--figpath', required=True, help='Path to where figures will be stored (end with /)')
    parser.add_argument('-n', '--countspath', required=True, help='Path to location where file that stores cDNA, gDNA, and normalized counts are stored.')
    parser.add_argument('-g', '--gdna', required=True, help='Gdna file with full path')
    parser.add_argument('-c', '--cdna', required=True, help='Cdna file with full path')
    parser.add_argument('-u', '--umi', action='store_true', help='Whether an umi is present in fastq at end of header')
    args = vars(parser.parse_args())

    FIG_PATH = args['figpath']
    COUNTS_PATH = args['countspath']
    cdna_reads = args['cdna']
    gdna_reads = args['gdna']
    umi_present = args['umi']

    FIG_PRE = '_'.join(cdna_reads.split('.')[0].split('_')[1:])
     
    # TODO turn collect reads into a separate file and use a snakemake file to run this
    # Generate filtered paired end fastqs
    #collect_reads([tag_file], [TAG_MARKER], [10], [tag_out])
    #collect_reads([tag_file, umi_file], [TAG_MARKER, UMI_MARKER], [10, 7], [tag_out, umi_out])

    # read in each umi_tag from header. After the initial first run of stats
    # have the choice to look at whether we want to discover why there were
    # error in some of the reverse reads by looking at UMI replicates
    tags = {}
    for line in open(TAGS, 'r'): 
        tags[line.rstrip('\n')] = True

    assert len(tags.keys()) == 7500

    cdna = grab_umi_tag(open(cdna_reads, 'r'), tags)
    gdna = grab_umi_tag(open(gdna_reads, 'r'), tags)

    # STAT 1 total reads post filtering, gDNA is most important stat here
    print("Post Filter Reads:")
    print("cDNA:", cdna.shape[0])
    print("gDNA:", gdna.shape[0])

    # STAT 2 total tags captured in gDNA
    print("Total unique tags:")
    print("cDNA:", cdna['Tag'].nunique())
    print("gDNA:", gdna['Tag'].nunique())

    # STAT 3 cDNA and gDNA counts post UMI collapse
    if umi_present:
        print("cDNA", end="\t")
        collapse_cdna = collapse_umis(cdna).reset_index()
    collapse_cdna = cdna.groupby(["Tag"]).size().reset_index()
    collapse_cdna.columns = ["Tag", "cDNA_counts"]

    if umi_present:
        print("gDNA", end="\t")
        collapse_gdna = collapse_umis(gdna).reset_index()
    collapse_gdna = gdna.groupby(["Tag"]).size().reset_index()
    collapse_gdna.columns = ["Tag", "gDNA_counts"]

    # Merge dataframes and create normalized counts
    cdna_gdna = collapse_gdna.merge(collapse_cdna, on='Tag', how='outer')
    cdna_gdna.fillna(0, inplace=True)
    print("Shared cDNA gDNA tags:", cdna_gdna[(cdna_gdna['cDNA_counts'] > 0) & (cdna_gdna['gDNA_counts'] > 0)].shape[0])
    cdna_gdna['cDNA/gDNA'] = cdna_gdna.apply(lambda row: row['cDNA_counts']/row['gDNA_counts'] if row['gDNA_counts'] > 0 else 0, axis=1)

    # Filter out low gDNA counts
    # TODO figure out what a good number is for noise removal?
    print(cdna_gdna.shape)
    filtered_cdna_gdna = cdna_gdna[cdna_gdna['gDNA_counts'] >= 100]
    print(filtered_cdna_gdna.shape)

    # output graphs that compare each tag and their respective counts
    normalize_tag_comparison(filtered_cdna_gdna, FIG_PATH, FIG_PRE, True)
    
    # create variance distributions of each original and random population
    true_var = variance(filtered_cdna_gdna) 
    rand_var = variance(filtered_cdna_gdna, rand=True)
  
    # used to zoom in on 0 and 0.5 for bins to prevent it from being like 300) 
    bins = np.linspace(0, 0.5, 50)
    bins = np.append(bins, np.linspace(0.5, 2.5))
    # plot histogram of variances
    fig, ax = plt.subplots(figsize=(7,7))
    ax.hist(true_var, bins=bins, alpha=0.5, label='True', color='b')
    ax.hist(rand_var, bins=bins, alpha=0.5, label='Random', color='r')
    ax.set_xlabel("Variance of Tag Replicates")
    ax.set_ylabel("Counts")
    ax.set_title("Distribution of True and Random Variances")
    ax.legend()
    plt.show()
    fig.tight_layout()
    fig.savefig(FIG_PATH + "variance_hist_" + FIG_PRE + ".png", bbox_inches='tight')

    # plot histogram of tags
    gdna_tag_counts = FIG_PATH + "counts_hist_" + gdna_reads.split('/')[-1].split('.')[0] + ".png"
    if not os.path.isfile(gdna_tag_counts):
        hist = collapse_gdna.hist(column="gDNA_counts", grid=False, bins=200)
        fig = hist[0][0].get_figure()
        fig.savefig(gdna_tag_counts)

    cdna_tag_counts = FIG_PATH + "counts_hist_" + cdna_reads.split('/')[-1].split('.')[0] + ".png"
    if not os.path.isfile(cdna_tag_counts):
        hist = collapse_cdna.hist(column="cDNA_counts", grid=False, bins=200)
        fig = hist[0][0].get_figure()
        fig.savefig(cdna_tag_counts)

    # write out cdna and gdna counts to create file
    cdna_gdna_file = COUNTS_PATH + "cdna_gdna_counts_" + FIG_PRE + ".txt"
    if not os.path.isfile(cdna_gdna_file):
        cdna_gdna.to_csv(cdna_gdna_file, sep='\t', index=False)
    

if __name__ == '__main__':
    main()
