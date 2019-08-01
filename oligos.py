import argparse
import sys
import random
import operator
import vcf
from itertools import imap

# Repeat all oligos 3 times each STR variant will have 5 permutations
#5'-ACTGGCCGCTTCACTG-var-GGTACCTCTAGA-tag-AGATCGGAAGAGCGTCG-3'
#5'-F1-var-KpnI-filler_seq-XbaI-tag-R1-3'   total size = 230 nt

KpnI = "GGTACC"
XbaI = "TCTAGA"

F1 = "ACTGGCCGCTTCACTG"
R1 = "AGATCGGAAGAGCGTCG"


# find hamming distasnce between two sequences
def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    ne = operator.ne
    return sum(imap(ne, seq1, seq2))


# Create initial list of eSTRs and indicate what strand they are present on
# eSTR in format chrom, start, end, motif
def gen_eSTRs(eSTRs, strands):
    strand_eSTRs = []
    for index, eSTR in enumerate(eSTRs):
        for strand in strands:
            if eSTR[4] == strand[1]:
                # check for forward strand or minus strand motif
                if strand[2] == '+':
                    strand_eSTRs.append([eSTR[0], eSTR[1], eSTR[2], eSTR[10]])
                else:
                    strand_eSTRs.append([eSTR[0], eSTR[1], eSTR[2], eSTR[11]])
    return strand_eSTRs


# 30 poly T
# 30 poly AC
# 20 random tetranucleotides
def gen_neg_cntrls(all_STRs, all_eSTRs):
    tetranucleotides = set()
    poly_AC = set()
    poly_T = set()
    negative_controls = []

    # Grab random tetranucleotide STRs for negative control
    for STR in all_STRs:
        if len(list(tetranucleotides)) < 20:
            if len(STR[3]) == 4:
                neg_cntrl = check_str(STR, all_eSTRs)
                if neg_cntrl:
                    tetranucleotides.add(neg_cntrl)
        if len(list(poly_AC)) < 30:
            if STR[3] == 'AC':
                neg_cntrl = check_str(STR, all_eSTRs)
                if neg_cntrl:
                    poly_AC.add(neg_cntrl)
        if len(list(poly_T)) < 30:
            if STR[3] == 'A':
                neg_cntrl = check_str(STR, all_eSTRs)
                if neg_cntrl:
                    poly_T.add(len(neg_cntrl)*'T')

    negative_controls.extend(list(tetranucleotides))
    negative_controls.extend(list(poly_AC))
    negative_controls.extend(list(poly_T))
    return negative_controls


# Check to see if STR is a currently existing eSTR
def check_str(STR, eSTRs):
    existing_eSTR = False
    for eSTR in eSTRs:
        if eSTR[0] == STR[0] and eSTR[1] == STR[1] and eSTR[2] == STR[2]:
            existing_eSTR = True
            break

    if not existing_eSTR:
        return ((int(STR[2])-int(STR[1]))/len(STR[3]))*STR[3]
    else:
        return
    

#TODO
def gen_eSTR_permutations(eSTRs, min_max_vcf):
    permutations = []

    # go over all eSTRs and determine all permutations
    for eSTR in eSTRs:
        for record in min_max_vcf.fetch(chrom=eSTR[0][3:], start=int(eSTR[1]), end=int(eSTR[2])):
            min_ref_del = 0
            max_ref_ins = 0

            for sample in record.samples:
                if sample.data[1] == None: continue
                GB = sample.data[1].split('|')
                for bp_dif in GB:
                    if int(bp_dif) > max_ref_ins:
                        max_ref_ins = int(bp_dif)
                    if int(bp_dif) < min_ref_del:
                        min_ref_del = int(bp_dif)

            # TODO
            #start = 
            #end = 

            # Determine whether we should use GTEX or vcf for start and end (record.start and record.end for vcf)
            min_repeat = (start-end+min_ref_del)//len(eSTR[3])
            max_repeat = (start-end+max_ref_ins)//len(eSTR[3])

            # check if there are no alternate alleles
            if not min_repeat == max_repeat:
                permutations.extend(['', min_repeat*eSTR[3], max_repeat*eSTR[3]])
            else:
                permutations.extend(['', min_repeat*eSTR[3], (min_repeat-1)*eSTR[3], (min_repeat+1)*eSTR[3], (min_repeat+2)*eSTR[3]])

            # one permutation before and after 
            if max_repeat - min_repeat == 1:
                permutations.extend([(min_repeat-1)*eSTR[3], (max_repeat+1)*eSTR[3]])
            # one permutation in between, one before
            elif max_repeat - min_repeat == 2:
                permutations.extend([(min_repeat-1)*eSTR[3], (min_repeat+1)*eSTR[3]])
            # 2 permuations in between try to make even distance
            else:
                permutatons.append((min_repeat + ((max_repeat - min_repeat)//3))*eSTR[3])
                permutations.append((max_repeat - ((max_repeat - min_repeat)//3))*eSTR[3])

    return permutations


# TODO do with Melissa, Alon, Catherine, and Sharona
def gen_fun():
    return


# Generate all oligos of eSTRs, negative controls, and the fun categories
# each oligo will contain three copies
def create_oligos(tags, filler, eSTRs=[], neg_controls=[], fun=[]):
    oligos = []
    var = []

    # Mark all types of var region sequences
    for eSTR in eSTRs:
        var.append((eSTR, 'eSTR'))
    for neg in neg_controls:
        var.append((neg, 'Negative Control'))
    for fun_str in fun:
        var.append((fun, 'For Fun'))

    # create three copies of each 230 nt oligo, every oligo has a unique tag
    t_ind = 0
    for seq in var:
        filler_len = 230 - (len(seq[0]) + len(tags[t_ind]) + 45)
        if filler_len > 0:
            filler_seq = filler[:filler_len]
        else:
            filler_seq = ''

        for i in range(3):
            oligos.append((seq[-1], F1 + seq[0] + KpnI + filler_seq + XbaI + tags[t_ind] + R1))
            t_ind += 1
        
    return oligos


# filter all currently existing oligos based on more than expected restriction enzyme sites
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
def output_oligos(oligos, output_file):
    with open(output_file, 'a') as output:
        output.write("Oligo_Type\tOligo\n")
        for oligo in oligos:
           output.write("%s\t%s\n"%(oligo[0], oligo[1]))
    return


def main():
    # All input files
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

    # file to output generated oligonucleotides too
    output_file = sys.argv[1]

    # Generate tags for each oligo nucleotide to be created
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
    #strand_eSTRs = gen_eSTRs(eSTRs[:450], strands)
    #str_permutations = gen_eSTR_permutations(strand_eSTRs, vcf_reader)
    negative_controls = gen_neg_cntrls(all_strs, eSTRs)
    #fun = gen_fun()

    oligos = create_oligos(all_tags, filler, eSTRs=str_permutations, neg_controls=negative_controls)
    filtered_oligos = filter_oligos(oligos)
    output_oligos(filtered_oligos, output_file)


if __name__ == '__main__':
    main()


