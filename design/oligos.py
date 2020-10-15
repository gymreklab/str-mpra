import sys
import random
import operator
import vcf
import pysam
import re
import math
import pandas as pd

debug_stretch = False

# Repeat all oligos 3 times each STR variant will have 5 permutations
#5'-F1-var-KpnI-filler_seq-XbaI-tag-R1-3'   total size = 230 nt

# Restriction Enzymes 
KpnI = "GGTACC"
XbaI = "TCTAGA"

# Primers 
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
# if the two sequences have a hamming distance less than the threshold
# return True otherwise return False
def hamming_distance(seq1, seq2, thresh):
    assert len(seq1) == len(seq2)
    dist = 0
    too_close = True
    for ind in range(len(seq1)):
        if seq1[ind] != seq2[ind]:
            dist += 1
        if dist >= thresh:
            too_close = False
            break
    return too_close


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


# Generate negative controls that are STRs not present in the eSTR GTEX data based on motifs inputted 
# /storage/resources/dbase/human/hg19/hg19.hipstr_reference_withmotif_stranded.bed IMPORTANT, helps check strand 
def negative_controls(all_STRs, all_eSTRs, controls):
    neg_cntrls = []
    strand_data = pd.read_csv('/storage/resources/dbase/human/hg19/hg19.hipstr_reference_withmotif_stranded.bed', 
                              names=['chrom', 'start', 'end', 'period', '+', '-'], sep='\t') 

    # Get current existing list of Motifs for negative control
    for STR in all_STRs:
        if all([control["total"] == control["current"] for control in controls]): break
        cur_strand = strand_data[(strand_data['chrom'] == STR[0]) & (strand_data['start'] == int(STR[1]))]['+'].item()

        for control in controls:
            # Already found all negative controls for current repeat unit
            if control["total"] == control["current"]: continue
            neg_cntrl = {}
            # search for control motif
            if re.search('^'+control["motif"]+'$', STR[3]):
                neg_cntrl = check_str(STR, all_eSTRs)
                if neg_cntrl:
                    # check if the Motif is on the + strand, if not then on - strand
                    if STR[3] == cur_strand:
                        neg_cntrl["strand"] = '+'
                    else:
                        neg_cntrl["strand"] = '-'

            # search for control reverse complement
            elif (not '[' in control["motif"]) and re.search('^'+reverse_complement(control["motif"])+'$', STR[3]):
                neg_cntrl = check_str(STR, all_eSTRs)
                if neg_cntrl:
                    neg_cntrl["motif"] = control["motif"]
                    # check reverse complement motif if matches then motif is on - strand
                    if STR[3] == cur_strand:
                        neg_cntrl["strand"] = '-'
                    else:
                        neg_cntrl["strand"] = '+'

            # found a negative control add it to the final list
            if neg_cntrl:
                control["current"] += 1
                neg_cntrls.append(neg_cntrl)

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

            # update min and max alleles to be factored by motif length
            if min_ref < 0:
                min_ref = len(STR["motif"]) * math.ceil(min_ref*1.0/len(STR["motif"]))
            else:
                min_ref = len(STR["motif"]) * math.floor(min_ref*1.0/len(STR["motif"]))
            if max_ref < 0:
                max_ref = len(STR["motif"]) * math.ceil(max_ref*1.0/len(STR["motif"]))
            else:
                max_ref = len(STR["motif"]) * math.floor(max_ref*1.0/len(STR["motif"]))

            # Start and end of repeat section
            start = record.INFO['START']
            end = record.INFO['END']
            zero_allele = -1*(end-start+1)

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
                split = (max_ref - min_ref + 1)/3.0
                from_min = len(STR["motif"]) * math.floor((min_ref + split)/len(STR["motif"]))
                from_max = len(STR["motif"]) * math.ceil((max_ref - split)/len(STR["motif"]))
                repeat_diffs.extend([from_min, from_max])

            # Generate all alleles and check if we want context to fill whole 175 or not
            for repeat_diff in repeat_diffs:
                if variable_context:
                    seq = _create_seq(STR, ref_genome, start, end, repeat_diff, repeat_diff)

                    # generate where the sequence starts in the reference for the allele
                    max_repeat_length = (end-start+1)+repeat_diff
                    context_len = (175-max_repeat_length)/2.0
                    seq_start = start-(math.floor(context_len))
                else:
                    seq = _create_seq(STR, ref_genome, start, end, repeat_diff, max_ref)
                    
                    # generate where seq starts in the reference for the allele
                    max_repeat_length = (end-start+1)+max_ref
                    context_len = (175-max_repeat_length)/2.0
                    seq_start = start-(math.floor(context_len))

                if seq is not None:
                    alleles.append(_create_allele(STR, seq, repeat_diff, seq_start))
                else:
                    print("WARNING: NO sequence created")
            break
    return alleles


# Generate all alleles with the given repeat number, header, and oligo information
def _create_allele(STR, seq, repeat_number, seq_start):
    return {"chrom":STR.get("chrom", ''), "start_pos":seq_start, 
            "str_pos":STR.get("start", ''), "str_end":STR.get("end", ''), 
            "gene":STR.get("gene", ''), "motif": STR["motif"], "num_repeats":repeat_number, "seq":seq}


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
        repeat_diff = -1 * (((-1*repeat_diff)//len(motif)) * len(motif))
        return left_region+rep_region[-1*repeat_diff:]+right_region
    else:
        repeat_diff = (repeat_diff//len(motif)) * len(motif)
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
    return final_seq


# generate fun sequences which are variable genomic context, random sequence replacing STR, and replace motif
def gen_fun(fun_file):
    fun_data = pd.read_csv(fun_file, sep='\t')
    fun_data.columns = ["chrom", "str_pos", "num_repeats", "seq"]
    fun = fun_data.to_dict(orient='records')
    return fun


# Generate all oligos of eSTRs, negative controls, and the fun categories
# Flag is used for the label of each oligo generated to give meta information
# each oligo has the form Forward primer + repeat sequeunce and genomic context + restriction enzyme 1 + filler seq + restriction enzyme 2 + tag + reverse primer
def create_oligos(tags, filler, all_alleles, flags=[]):
    oligos = []
    assert len(all_alleles) == len(flags)

    # create three copies of each 230 nt oligo, every oligo has a unique tag
    t_ind = 0
    for var, flag in zip(all_alleles, flags):
        for seq in var:
            filler_len = 230 - (len(seq["seq"]) + len(tags[t_ind]) + 45)
            if filler_len > 0:
                filler_seq = filler[:filler_len]
            else:
                filler_seq = ''

            # marker to label seq when printing, use flag to indicate if marker is wanted
            if flag:
                marker = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(flag, seq["chrom"], str(seq.get("start_pos", '')), 
                                                 str(seq["str_pos"]), str(seq.get("str_end", '')), seq.get("gene", ''), 
                                                 seq.get("motif", ''), str(seq["num_repeats"]))

            for i in range(3):
                oligos.append({"label":marker, "seq":' '.join([F1, seq["seq"], KpnI, filler_seq, XbaI, tags[t_ind], R1])})
                t_ind += 1
        
    return oligos


# filter all currently existing oligos based on more than expected restriction enzyme sites
def filter_oligos(oligos, restriction_sites):
    filtered_oligos = []

    # iterate over all oligo and check if exists more restriction sites than expected
    for oligo in oligos:
        sequence = oligo["seq"].replace(' ', '')
        assert len(sequence) == 230
        filter_seq = False
        extra_site_locs = []
        for site, threshold in restriction_sites:
            if len(re.findall(site, sequence)) > threshold:
                filter_seq = True
                extra_site_locs.extend(re.findall(site, sequence))
        if not filter_seq:
            filtered_oligos.append(oligo)
        else:
            print("Filtered_oligo %s %s"%(oligo["label"], oligo["seq"]))
            print("Extra sites:", extra_site_locs)

    return filtered_oligos


# Write all filtered oligos to a file
def output_oligos(oligos, output_file):
    with open(output_file, 'a') as output:
        output.write("Oligo_Type\tChrom\tStart_pos\tSTR_pos\tSTR_end\tGene\tMotif\tRepeat_Number\tFoward_Primer_1 Variant_Sequence Restriction_Enzyme_1 Filler_Sequence Restriction_Enzyme_2 Tag Reverse_Primer_1\n")
        for oligo in oligos:
           output.write("%s\t%s\n"%(oligo["label"], oligo["seq"]))
    return


def main():
    # Input Data 
    ref_genome = pysam.Fastafile('/storage/resources/dbase/human/hg19/hg19.fa')
    tag_file = open('./permutation_tags.txt', 'r')
    fun_file = './fun_loci/STRMPRA_FunLoci.txt'
    vcf_reader = vcf.Reader(open('/storage/szfeupe/Runs/650GTEx_estr/Filter_Merged_STRs_All_Samples_New.vcf.gz', 'rb'))
    strand_file = '/storage/mlamkin/projects/eSTR-data/gencode.v7.tab'

    eSTR_file_path = '/storage/mlamkin/projects/eSTR-data/eSTRGtex_DatasetS1.csv'
    all_strs_file = '/storage/mlamkin/projects/eSTR-data/all_analyzed_strs_v2.tab'
    filler = open('./restriction_filler.txt', 'r').readline().rstrip('\n')

    # input data
    all_strs = []
    eSTRs = []
    strands = []
    all_tags = []
    oligos = []

    # file to output generated oligonucleotides too
    output_file = sys.argv[1]

    # Grab tags from pregenerated file for each oligo to be created
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
    strand_eSTRs = find_eSTRs(eSTRs[:430], strands)
    eSTR_alleles = create_alleles(strand_eSTRs, vcf_reader, ref_genome)
    # pass in dictionary of STRs as regex to determine what controls are wanted
    neg_cntrls = negative_controls(all_strs, eSTRs, [{'motif':'T', 'total':30, 'current':0}, {'motif':'AC','total':30, 'current':0},{'motif': '[ACGT]{4}','total':20, 'current':0}])
    neg_alleles = create_alleles(neg_cntrls, vcf_reader, ref_genome)
    fun_alleles = gen_fun(fun_file)

    # Generate Oligos with all alleles found
    all_alleles = [neg_alleles, eSTR_alleles, fun_alleles]
    flags = ['Negative_Control', 'eSTR', 'Fun']
    oligos.extend(create_oligos(all_tags, filler, all_alleles, flags=flags))

    # restriction sites expected to be in sequence
    restriction_sites = [(KpnI, 1), (XbaI, 1), (r'GGCC[ACGT]{5,5}GGCC', 0)]

    # filter and output oligos
    filtered_oligos = filter_oligos(oligos, restriction_sites)
    output_oligos(filtered_oligos, output_file)


if __name__ == '__main__':
    main()

