import argparse
import sys
import random
import operator
import vcf
import pysam
import math
from itertools import imap

# Repeat all oligos 3 times each STR variant will have 5 permutations
#5'-ACTGGCCGCTTCACTG-var-GGTACCTCTAGA-tag-AGATCGGAAGAGCGTCG-3'
#5'-F1-var-KpnI-filler_seq-XbaI-tag-R1-3'   total size = 230 nt

#TODO Make the Restriction enzymes and Primers a separate file to load in by the user
KpnI = "GGTACC"
XbaI = "TCTAGA"

F1 = "ACTGGCCGCTTCACTG"
R1 = "AGATCGGAAGAGCGTCG"

# find reverse complement of a sequence
def reverse_complement(seq):
    comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    rev_comp = ''
    for nuc in seq:
        rev_comp = rev_comp + comp[nuc]
    return rev_comp[::-1]


# find hamming distasnce between two sequences
def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    ne = operator.ne
    return sum(imap(ne, seq1, seq2))


# Create initial list of eSTRs and indicate what strand they are present on
# eSTR in format chrom, start, end, motif
def find_eSTRs(eSTRs, strands):
    strand_eSTRs = []
    for index, eSTR in enumerate(eSTRs):
        for strand in strands:
            if eSTR[4] == strand[1]:
                # check for forward strand or minus strand motif
                if strand[2] == '+':
                    strand_eSTRs.append({"chrom":eSTR[0], "start":eSTR[1], 
                                         "end":eSTR[2], "gene":eSTR[4], "motif":eSTR[10]})
                else:
                    strand_eSTRs.append({"chrom":eSTR[0], "start":eSTR[1], 
                                         "end":eSTR[2], "gene":eSTR[4], "motif":eSTR[11]})
    return strand_eSTRs


# 30 poly T
# 30 poly AC
# 20 random tetranucleotides
# TODO Make more user friendly so you can change what kind of negative controls you want and how many
def negative_controls(all_STRs, all_eSTRs):
    tetranucleotides = []
    poly_AC = []
    poly_T = []
    neg_cntrls = []

    # Grab random tetranucleotide STRs for negative control
    for STR in all_STRs:
        if len(tetranucleotides) < 20:
            if len(STR[3]) == 4:
                neg_cntrl = check_str(STR, all_eSTRs)
                if neg_cntrl:
                    tetranucleotides.append(neg_cntrl)
        if len(poly_AC) < 30:
            if STR[3] == 'AC':
                neg_cntrl = check_str(STR, all_eSTRs)
                if neg_cntrl:
                    poly_AC.append(neg_cntrl)
        if len(poly_T) < 30:
            if STR[3] == 'A':
                # only has poly A's and no way to check strand so alter to T
                neg_cntrl = check_str(STR, all_eSTRs)
                if neg_cntrl:
                    neg_cntrl["motif"] = 'T'
                    poly_T.append(neg_cntrl)

    neg_cntrls.extend(tetranucleotides)
    neg_cntrls.extend(poly_AC)
    neg_cntrls.extend(poly_T)
    return neg_cntrls


# Check to see if STR is a currently existing eSTR
def check_str(STR, eSTRs):
    existing_eSTR = False
    for eSTR in eSTRs:
        if eSTR[0] == STR[0] and eSTR[1] == STR[1] and eSTR[2] == STR[2]:
            existing_eSTR = True
            break

    if not existing_eSTR:
        return {"chrom":STR[0], "start":STR[1], "end":STR[2], "motif":STR[3]}
    else:
        return
    

#TODO ADD HEADER INTO ALLELE
# THINK ABOUT CHANGING THE LISTS TO DICTIONARIES FOR CLEANER CODE INSTEAD OF HAVING TO FORCE THE MOTIF TO BE IN THE END
# ONLY NEED TO HAVE CHROM START END MOTIF BUT CAN CONTAIN OTHER WHICH WE CAN ACCOMODATE FOR HEADER IN OLIGO
def create_alleles(STRs, min_max_vcf, ref_genome):
    alleles = []

    # grab gene name from eSTRs for marker for each STR (if neg control wont have gene name)
    # go over all eSTRs and determine all permutations
    for STR in STRs:
        for record in min_max_vcf.fetch(chrom=STR["chrom"][3:], start=int(STR["start"]), end=int(STR["end"])):
            # min and max bp differences for all samples vs reference to get range
            min_ref = 0
            max_ref = 0

            # iterate through each sample and grab each bp difference from reference
            for sample in record.samples:
                if sample.data[1] == None: continue
                GB = [int(sample_allele) for sample_allele in sample.data[1].split('|')]
                if min(GB) < min_ref: min_ref = min(GB)
                if max(GB) > max_ref: max_ref = max(GB)

            # Start and end of repeat section
            start = record.INFO['START']
            end = record.INFO['END']

            # Find the min and max alleles of each STR
            min_repeat = (end-start+min_ref)//len(STR["motif"])
            max_repeat = (end-start+max_ref)//len(STR["motif"])

            # length of genomic context to bring total length of max allele STR to 175 nucleotides
            context_len = (175-(max_repeat*len(STR["motif"])))/2.0

            # Make sure context is from the proper strand TODO (Maybe store strand)
            l_context = ref_genome.fetch(region='%s:%i-%i'%(STR["chrom"], start-(math.floor(context_len)), start-1)).upper()
            r_context = ref_genome.fetch(region='%s:%i-%i'%(STR["chrom"], end+1, end+math.ceil(context_len))).upper()

            if STR.get("strand", None) == '-':
                r_context = reverse_complement(l_context)
                l_context = reverse_complement(r_context)

            # check if there are no alternate alleles (Throw flag because this should not happen)
            if not min_repeat == max_repeat:
                alleles.extend([_create_allele(STR, l_context + repeat*STR["motif"] + r_context, repeat) for repeat in [0, min_repeat, max_repeat]])
            else:
                print("sample found with no allelic differences: %s %i %i"%(STR["chrom"], start, end))

            # one permutation before and after 
            if max_repeat - min_repeat == 1:
                alleles.extend([_create_allele(STR, l_context + repeat*STR["motif"] + r_context, repeat) for repeat in [min_repeat-1, max_repeat+1]])
            # one allele in between, one before
            elif max_repeat - min_repeat == 2:
                alleles.extend([_create_allele(STR, l_context + repeat*STR["motif"] + r_context, repeat) for repeat in [min_repeat-1, min_repeat+1]])
            # 2 alleles in between try to make even distance
            else:
                alleles.extend([_create_allele(STR, l_context + repeat*STR["motif"] + r_context, repeat) \
                               for repeat in [(min_repeat + ((max_repeat - min_repeat)//3)), (max_repeat - ((max_repeat - min_repeat)//3))]])
            
    return alleles


# Generate all alleles with the given repeat number, header, and oligo information 
def _create_allele(STR, oligo, repeat_number):
    return {"chrom":STR.get("chrom", ''), "pos":STR.get("start", ''),
            "gene":STR.get("gene", ''), "num_repeats":repeat_number,
            "oligo":oligo}


# TODO do with Melissa, Alon, Catherine, and Sharona
def gen_fun():
    
    return


# Generate all oligos of eSTRs, negative controls, and the fun categories
# Make function more generic, just do it with data given
# Header is made up of everything except last col which contains STR allele
# Can have flags to indicate whether it is an eSTR, negative control, or fun oligo
# TODO UPDATE WItH NEW HEADER FORMAT IN DICTIONARY
def create_oligos(tags, filler, var, flag=''):
    oligos = []

    # create three copies of each 230 nt oligo, every oligo has a unique tag
    t_ind = 0
    for seq in var:
        filler_len = 230 - (len(seq[-1]) + len(tags[t_ind]) + 45)
        if filler_len > 0:
            filler_seq = filler[:filler_len]
        else:
            filler_seq = ''

        # marker to label seq when printing, use flag to indicate if marker is wanted
        if flag:
            marker = flag + '_' + '_'.join([])

        for i in range(3):
            oligos.append((marker, F1 + seq[-1] + KpnI + filler_seq + XbaI + tags[t_ind] + R1))
            t_ind += 1
        
    return oligos


# filter all currently existing oligos based on more than expected restriction enzyme sites
# TODO Add in user defined restriction enzymes and reverse complement of said enzymes
def filter_oligos(oligos):
    filtered_oligos = []
    comp_motifs = [KpnI, XbaI]

    # iterate over all oligo and check if more than 2 potential cut sites
    # if more than 2 than filter that oligo out
    for oligo in oligos:
        thresh = 0
        matching = []
        for motif in comp_motifs:
            assert len(motif) <= len(oligo[1])
            for i in range(len(oligo[1])-len(motif)+1):
                if hamming_distance(motif, oligo[1][i:i+len(motif)]) < 2:
                    matching.append(oligo[1][i:i+len(motif)])
                    thresh += 1
        if thresh == 2:
            filtered_oligos.append(oligo)
        else:
            print(oligo, matching)

    return filtered_oligos


# Write all filtered oligos to a file
# TODO UPDATE WITH NEW FORMAT
def output_oligos(oligos, output_file):
    with open(output_file, 'a') as output:
        output.write("Oligo_Type\tChrom\tPos\tGene\tRepeatNumber\tOligo\n")
        for oligo in oligos:
           output.write("%s\t%s\n"%(oligo[0], oligo[1]))
    return


def main():
    # All input files
    ref_genome = pysam.Fastafile('/storage/resources/dbase/human/hg19/hg19.fa') # TODO make this a user input(hg19 fasta path)
    tag_file = open('/storage/mlamkin/projects/str-mpra/permutation_tags.txt', 'r')
    vcf_reader = vcf.Reader(open('/storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz', 'r'))
    strand_file = '/storage/mlamkin/projects/eSTR-data/gencode.v7.tab'
    eSTR_file_path = '/storage/mlamkin/projects/eSTR-data/eSTRGtex_DatasetS1.csv'
    all_strs_file = '/storage/mlamkin/projects/eSTR-data/all_analyzed_strs_v2.tab'
    filler = open('/storage/mlamkin/projects/str-mpra/restriction_filler.txt', 'r').readline().rstrip('\n')

    # input data
    all_strs = []
    eSTRs = []
    strands = []
    all_tags = []
    oligos = []

    # file to output generated oligonucleotides too
    output_file = sys.argv[1]

    # Generate tags for each oligo nucleotide to be created
    # TODO make user specified amount of tags to be taken
    for line in tag_file:
        all_tags.append(line.rstrip('\n'))

    # read in GTEX eSTRs
    with open(eSTR_file_path) as eSTR_file:
        eSTR_file.readline()
        for line in eSTR_file:
            eSTRs.append(line.rstrip('\n').split(','))

    # Read in all STRs with their strand and gene association
    with open(strand_file) as strand:
        next(strand)
        for line in strand:
            strands.append(line.rstrip('\n').split('\t'))

    # Read in all strs available foe negative control
    for line in open(all_strs_file):
        all_strs.append(line.rstrip('\n').split('\t')) 

    # create list of eSTRs
    # TODO maybe make list of eSTRs to take a user specified option
    strand_eSTRs = find_eSTRs(eSTRs[:1], strands)
    eSTR_alleles = create_alleles(strand_eSTRs, vcf_reader, ref_genome)
    # TODO ENSURE THAT STRAND REVERSE COMPLEMENT WORKS
    neg_cntrls = negative_controls(all_strs, eSTRs)
    neg_alleles = create_alleles(neg_cntrls, vcf_reader, ref_genome)
    #fun = gen_fun() TODO GENERATE FUN SEQUENCES

    # TODO Should make sure you dont have to manually input these but they will appear based on what user wants 
    for seqs in [(neg_alleles, 'Negative Control'), (eSTR_alleles, 'eSTR'), (fun_alleles, 'Fun')]:
        oligos.extend(create_oligos(all_tags, filler, seqs[0], flag=seqs[1]))

    filtered_oligos = filter_oligos(oligos)
    output_oligos(filtered_oligos, output_file)


if __name__ == '__main__':
    main()


