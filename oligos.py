import argparse
import sys
import random
import operator
import vcf
import pysam
import re
import math

debug_stretch = False

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
                    strand_eSTRs.append({"chrom":eSTR[0], "start":eSTR[1], "strand":strand[2],
                                         "end":eSTR[2], "gene":eSTR[4], "motif":eSTR[10]})
                else:
                    strand_eSTRs.append({"chrom":eSTR[0], "start":eSTR[1], "strand":strand[2],
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
    

# Generate 4 alleles per STR each of different length and one with only genomic context surrounding the STR in the reference
def create_alleles(STRs, min_max_vcf, ref_genome, variable_context=False):
    alleles = []

    assert isinstance(STRs, list)
    # grab gene name from eSTRs for marker for each STR (if neg control wont have gene name)
    # go over all eSTRs and determine all permutations
    for STR in STRs:
        for record in min_max_vcf.fetch(chrom=STR["chrom"][3:], start=int(STR["start"])):
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
            zero_allele = -1*(end-start+1) # TODO check this....

            # Hold repeat lengths for alleles as diff from ref
            repeat_diffs = []

            # check if there are no alternate alleles (Throw flag because this should not happen)
            if not min_ref == max_ref:
                repeat_diffs.extend([zero_allele, min_ref, max_ref])
            else:
                print("sample found with no allelic differences: %s %i %i"%(STR["chrom"], start, end))

            # two alleles before min_repeat number
            if max_ref - min_ref == len(STR["motif"]): 
                repeat_diffs.extend([min_ref-len(STR["motif"]), min_ref-2*len(STR["motif"])])
            # one allele in between, one before
            elif max_ref - min_ref == 2*len(STR["motif"]): 
                repeat_diffs.extend([min_ref-len(STR["motif"]), min_ref+len(STR["motif"])])
            # 2 alleles in between try to make even distance
            else: 
                repeat_diffs.extend([(min_ref + ((max_ref - min_ref)//3)), (max_ref - ((max_ref - min_ref)//3))])

            # Generate all alleles and check if we want context to fill whole 175 or not
            for repeat_diff in repeat_diffs:
                print("Repeat diff %s %s %s"%(repeat_diff, STR.get("gene","NA"), STR.get("strand","NA")))
                if variable_context:
                    seq = _create_seq(STR, ref_genome, start, end, repeat_diff, repeat_diff)
                else:
                    seq = _create_seq(STR, ref_genome, start, end, repeat_diff, max_ref)
                if seq is not None:
                    alleles.append(_create_allele(STR, seq, repeat_diff))
                else:
                    print("WARNING: NO sequence created")
            break
    return alleles


# Generate all alleles with the given repeat number, header, and oligo information 
def _create_allele(STR, seq, repeat_number):
    return {"chrom":STR.get("chrom", ''), "pos":STR.get("start", ''),
            "gene":STR.get("gene", ''), "num_repeats":repeat_number,
            "seq":seq}

def _is_perfect(sequence, motif):
    if debug_stretch: print("checking %s %s"%(sequence, motif))
    perf_seq = (motif*math.ceil(len(sequence)*1.0/len(motif)))[0:len(sequence)]
    return perf_seq==sequence

def _longest_repeat(sequence, motif):
    if debug_stretch: print("%s: %s"%(sequence, motif))
    rotations = []
    for i in range(len(motif)):
        rotations.append(motif[i:]+motif[0:i])
        
    longest_stretch_start = -1
    longest_stretch = -1
    best_rotation = None

    for rot in rotations:
        offset = 0
        loc = sequence.find(rot)
        while True:
            stretch_length = len(rot)
            while True:
                if offset+loc+stretch_length >= len(sequence): break
                if _is_perfect(sequence[offset+loc:offset+loc+stretch_length+1], rot):
                    stretch_length += 1
                else: break
            if stretch_length > longest_stretch:
                longest_stretch_start = loc+offset
                longest_stretch = stretch_length
                best_rotation = rot
            offset = offset+loc+stretch_length
            x = sequence[offset:].find(rot)
            if x == -1:
                break
            else: loc = x
            if debug_stretch: print("rot: %s loc: %s offset: %s x: %s"%(rot, loc, offset, x))
    if debug_stretch: print("Answer: %s %s %s"%(sequence[longest_stretch_start: longest_stretch_start+longest_stretch], longest_stretch_start, longest_stretch))
    return longest_stretch_start, longest_stretch, best_rotation

def _get_perfect_stretch(best_rotation, repeat_diff):
    motif_len = len(best_rotation)
    perf_seq = best_rotation*math.ceil(repeat_diff*1.0/motif_len)
    return perf_seq[0:repeat_diff]

def _get_repeat_seq(ref_repeat, repeat_diff, motif):
    start_offset, stretch_len, best_rotation = _longest_repeat(ref_repeat, motif)
    left_region = ref_repeat[0:start_offset]
    rep_region = ref_repeat[start_offset:start_offset+stretch_len]
    right_region = ref_repeat[start_offset+stretch_len:]

    # Special cases
    if start_offset == -1:
        print("WARNING: could not find %s in %s"%(motif, ref_repeat))
        return None
    if repeat_diff + len(ref_repeat) < 0:
        print("WARNING: repeat diff %s less than repeat length %s"%(repeat_diff, len(ref_repeat)))
        return None
    if repeat_diff + len(rep_region) < 0:
        print("WARNING: repeat diff %s less than repeat length %s. Deleting non-perfect repeat."%(repeat_diff, len(rep_region)))
        return ref_repeat[-1*repeat_diff:]

    # Add or delete from perfect stretches
    if repeat_diff == 0:
        return ref_repeat
    elif repeat_diff < 0:
        return left_region+rep_region[-1*repeat_diff:]+right_region
    else:
        new_chunk = _get_perfect_stretch(best_rotation, repeat_diff)
        return left_region+new_chunk+rep_region+right_region

def _create_seq(STR, ref_genome, start, end, repeat_diff, max_ref):
    max_repeat_length = (end-start+1)+max_ref
    context_len = (175-max_repeat_length)/2.0
    
    # genomic context left of the motif and right of the motif (reverse complement and switched if on opposite strand)
    l_context = ref_genome.fetch(region='%s:%i-%i'%(STR["chrom"], start-(math.floor(context_len)), start-1)).upper()
    r_context = ref_genome.fetch(region='%s:%i-%i'%(STR["chrom"], end+1, end+math.ceil(context_len))).upper()
    ref_repeat = ref_genome.fetch(region='%s:%i-%i'%(STR["chrom"], start, end)).upper()


    if STR.get("strand", None) == '-':
        newr_context = reverse_complement(l_context)
        l_context = reverse_complement(r_context)
        r_context = newr_context
        ref_repeat = reverse_complement(ref_repeat)

    repeat_seq = _get_repeat_seq(ref_repeat, repeat_diff, STR["motif"])
    if repeat_seq is None: return None
    
    final_seq = l_context + repeat_seq + r_context

    #### debugging
    print("########")
    refseq = ref_genome.fetch(region='%s:%i-%i'%(STR["chrom"], start-(math.floor(context_len)), end+math.ceil(context_len))).upper()
    if STR.get("strand", None) == '-':
        print(reverse_complement(refseq))
    else: print(refseq)
    print("%s: %s"%(final_seq, len(final_seq)))
    print("########")
    #### debugging

    return final_seq

# generate fun sequences which are variable genomic context, random sequence replacing STR, and replace motif
def gen_fun(vcf, ref_genome):
    # Dinucleotide regions do TG, AG, AAC * 4 for each so 12 
    # Random sequences (replace same stretch of AC but with random nucleotides) do this 8 times
    # normal region but include extra genomic context so there is no filler inbetween
    # FUN REGIONS:
    nucs = ['A', 'C', 'G', 'T']
    alt_allele_di = ['TG', 'AG', 'AAC']
    alt_allele_six = ['TGCGAG', 'CCGTCA', 'ACGC']
    RFT1 = {"chrom":"chr3", "start":'53128363', "motif":"AC"}
    APEH = {"chrom":"chr3", "start":'49711229', "motif":"ACGCTC"}
    DISP2 = {"chrom":"chr15", "start":'40643044', "motif":"AC"}
    GLYCTK = {"chrom":"chr3", "start":'52341945', "motif":"GT"}
    fun = []

    for locus in [RFT1, DISP2, GLYCTK, APEH]:
        for record in vcf.fetch(chrom=locus["chrom"][3:], start=int(locus["start"])):
            # Start and end of repeat section
            start = record.INFO['START']
            end = record.INFO['END']

            # repeat number in reference
            ref_repeat = ((end-start)//2)
            context_len = (175-(ref_repeat*len(locus["motif"])))/2.0

            l_context = ref_genome.fetch(region='%s:%i-%i'%(locus["chrom"], start-(math.floor(context_len)), start-1)).upper()
            r_context = ref_genome.fetch(region='%s:%i-%i'%(locus["chrom"], end+1, end+math.ceil(context_len))).upper()
             
            # determine what list of alternate alleles to use
            if len(locus["motif"]) == 2: alt_allele = alt_allele_di
            else: alt_allele = alt_allele_six

            # add 4 of each alternate allele
            for alt in alt_allele:
                alt_context_len = (175-(ref_repeat*len(alt)))/2.0
                alt_l_context = ref_genome.fetch(region='%s:%i-%i'%(locus["chrom"], start-(math.floor(alt_context_len)), start-1)).upper()
                alt_r_context = ref_genome.fetch(region='%s:%i-%i'%(locus["chrom"], end+1, end+math.ceil(alt_context_len))).upper()
                for i in range(4):
                    fun.append({"chrom":locus["chrom"], "pos":locus["start"], 
                                "num_repeats":ref_repeat,
                                "seq":alt_l_context+ref_repeat*alt+alt_r_context})
            # add 8 replacing the STR with random sequence holding the context constant
            for i in range(8):
                fun.append({"chrom":locus["chrom"], "pos":locus["start"], 
                            "num_repeats":"random_seq",
                            "seq":l_context+"".join([random.choice(nucs) for i in range(ref_repeat*len(locus["motif"]))])+r_context})
            break

    # last 4-5 have alternating repeats with different genomic context adding to 175
    fun.extend(create_alleles([DISP2, GLYCTK], vcf, ref_genome, variable_context=True)) 
    return fun


# Generate all oligos of eSTRs, negative controls, and the fun categories
# Flag is used for the label of each oligo generated to give meta information
# each oligo has the form Forward primer + repeat sequeunce and genomic context + restriction enzyme 1 + filler seq + restriction enzyme 2 + tag + reverse primer
def create_oligos(tags, filler, var, flag=''):
    oligos = []

    # create three copies of each 230 nt oligo, every oligo has a unique tag
    t_ind = 0
    for seq in var:
        filler_len = 230 - (len(seq["seq"]) + len(tags[t_ind]) + 45)
        if filler_len > 0:
            filler_seq = filler[:filler_len]
        else:
            filler_seq = ''

        # marker to label seq when printing, use flag to indicate if marker is wanted
        if flag:
            marker = "%s\t%s\t%s\t%s\t%s"%(flag, seq["chrom"], 
                        str(seq["pos"]), seq.get("gene", ''), str(seq["num_repeats"]))

        for i in range(3):
            oligos.append({"label":marker, "seq":' '.join([F1, seq["seq"], KpnI, filler_seq, XbaI, tags[t_ind], R1])})
            t_ind += 1
        
    return oligos


# filter all currently existing oligos based on more than expected restriction enzyme sites
# TODO Add in user defined restriction enzymes and reverse complement of said enzymes and sequences to scan with regex
# Label each enzyme with the number of times we expect it to be in the sequence
def filter_oligos(oligos):
    filtered_oligos = []

    # restriction sites expected to be in sequence only once
    restriction_sites = [(KpnI, 1), (XbaI, 1), (r'GGCC[ACGT]{5,5}GGCC', 0)]

    # iterate over all oligo and check if exists more restriction sites than expected
    for oligo in oligos:
        sequence = oligo["seq"].replace(' ', '')
        filter_seq = False
        for site, threshold in restriction_sites:
            if len(re.findall(site, sequence)) > threshold:
                filter_seq = True
        if not filter_seq:
            filtered_oligos.append(oligo)

    return filtered_oligos


# Write all filtered oligos to a file
def output_oligos(oligos, output_file):
    with open(output_file, 'a') as output:
        output.write("Oligo_Type\tChrom\tPos\tGene\tRepeatNumber\tFoward_Primer_1 Variant_Sequence Restriction_Enzyme_1 Filler_Sequence Restriction_Enzyme_2 Tag Reverse_Primer_1\n")
        for oligo in oligos:
           output.write("%s\t%s\n"%(oligo["label"], oligo["seq"]))
    return


def main():
    # TODO use argparse to set user defined paths for the majority of these values 
    # All input files
    ref_genome = pysam.Fastafile('/storage/resources/dbase/human/hg19/hg19.fa') # TODO make this a user input(hg19 fasta path)
    tag_file = open('/storage/mlamkin/projects/str-mpra/permutation_tags.txt', 'r')
    vcf_reader = vcf.Reader(open('/storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz', 'rb'))
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
    # TODO maybe make list of eSTRs to take a user specified option to determine how many wanted
    strand_eSTRs = find_eSTRs(eSTRs[:430], strands)
    eSTR_alleles = create_alleles(strand_eSTRs, vcf_reader, ref_genome)
    neg_cntrls = negative_controls(all_strs, eSTRs)
    neg_alleles = create_alleles(neg_cntrls, vcf_reader, ref_genome)
    fun_alleles = gen_fun(vcf_reader, ref_genome)
    # TODO Should make sure you dont have to manually input these but they will appear based on what user wants 
    for seqs in [(neg_alleles, 'Negative Control'), (eSTR_alleles, 'eSTR'), (fun_alleles, 'Fun')]:
        oligos.extend(create_oligos(all_tags, filler, seqs[0], flag=seqs[1]))

    filtered_oligos = filter_oligos(oligos)
    #output_oligos(filtered_oligos, output_file)


if __name__ == '__main__':
    main()

#seq="GTTTGTTTTGTTTGTTTGTTTGTTTGTTTGTTTGTTT"
#_longest_repeat(seq, "GTTT")

#seq="ATATATATATATATGTATAT"
#_longest_repeat(seq, "AT")
