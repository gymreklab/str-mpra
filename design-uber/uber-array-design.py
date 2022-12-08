#!/usr/bin/env python3
"""
Design uber STR array sequence for MPRA

Example command:
./uber-array-design.py \
   --reffa chr21.fa \
   --rptsbed chr21_test.bed \
   --out chr21_test


Questions:
- Check set of alternate motifs
- Filter redundant motifs. e.g those with length 0 will be same for all
"""


import argparse
import numpy as np
import os
import pyfaidx
import random
import sys

from trtools.utils import utils

# Global variables for oligo construction. These cannot change
FIVE_PRIME_ADAPT = 'ACTGGCCGCTTGACG'
GIBSON_ASISI = 'CACTGCGGCTCCTGCGATCGC'
GLOBAL_FILLER_SEQ = 'TAACCAGGCGTGTTAGCTGCTGTGCTGTCCTACGAGTAAACAGTAGAGTCCGTGGGCACGCGAGCACGGTGAGTCGACTCTGGCCTCATCACCATTTAGTTTGCGCAAGCGCTCTTTTTATAGGACCTGTCTTACATCCCTCATTAACGGAATCGATTACCGGCTAGCGTTGAAATGGAGAAACCGGCTTGCAGTCGAAA'
BSAI_RECOG = 'GGTCTC'
AsiSI_RECOG = 'GCGATCGC'
GIBSON_BSAI_CUT = 'TGTCGATCGCGTCGACGAAC'
PROBE_LEN = 230
MAX_STR_LEN = 80

# Other configuration variables that can be played with
TOTAL_LENGTHS = {1: 20, 2: 15, 3: 25, 4: 20, 5: 16, 6: 13}

# desired construct for each STR
# 5'adapter-(var + genomic context)-(gibson_asis)+(filler_bsai_recog)+(gibson_bsai_cut)

# Note: can change if we want to consider a different
# set of unique repeat units lengths.
# e.g. rather than 0, 1, 2, 3, etc. can do 0, 2, 4, 6, etc.
def GetAlleleLengths(num_rpt_lengths, repeat_unit_length):
	"""
	Determine the unique allele lengths to generate

	Arguments
	---------
	num_rpt_lengths: int
	   The number of unique lengths we want to generate
	repeat_unit_length : int
	   The length of the repeat unit

	Returns
	-------
	rpt_lengths : list of int
	   List of repeat lengths to consider
	"""
	return list(range(0, num_rpt_lengths))

# Susan: think about if we want to take the top 2?
# or randomly choose?
# Is this an appropriate set of alternate motifs?
def GetAlternateMotifs(exclude=None):
	"""
	Get a list of alternative motifs to consider
	We will return 2 each of di-, tri-, tetra-

	Arguments
	---------
	exclude : str (optional)
	   Optional list of motifs to exclude

	Returns
	-------
	motifs : list of str
	   List of alternate motifs to consider
	"""
	alt_di_motifs = ["AC", "AT", "AG", "GC"]
	alt_tri_motifs = ["CCG","AAT","ACG","AAG"]
	alt_tetra_motifs = ["AAAC","AGAT","AAAT","ACAT"]

	# Remove excluded motif from these lists upfront
	if exclude is not None:
		if exclude in alt_di_motifs: alt_di_motifs.remove(exclude)
		if exclude in alt_tri_motifs: alt_tri_motifs.remove(exclude)
		if exclude in alt_tetra_motifs: alt_tetra_motifs.remove(exclude)
	motifs = alt_di_motifs[0:2] + alt_tri_motifs[0:2] + alt_tetra_motifs[0:2]
	return motifs

def GenerateRandomSequence(seqlen, gcperc=0.5, maxtries=100):
	"""
	Generate random sequence. 
	Keep trying until we don't have a cut site

	Arguments
	---------
	seqlen : int
	  Length of the sequence to generate
	gcperc : float (optional)
	  GC percentage to target
	maxtries : int
	  Maximum times to try to generate a sequence 
	  that passes with no cut sites

	Returns
	-------
	random_seq : str
	  Random sequence

	Example
	-------
	GenerateRandomSequence(10, gcperc=0.0)
	> ATTATTATAT
	"""
	passed = False
	random_seq = ''
	numtries = 0
	while not passed:
		random_seq = ''.join(random.choices('ATCG', weights = ((1-gcperc)/2, (1-gcperc)/2, gcperc/2, gcperc/2), k=seqlen))
		if CheckCutSites(random_seq): passed = True
		numtries += 1
		if numtries > maxtries:
			sys.stderr.write("Failed to generate a random sequence with no cut sites. Aborting mission...\n")
			sys.exit(1)
	return random_seq
                   

def GetGC(seq):
	"""
	Compute the GC percentage of a sequence

	Arguments
	---------
	seq : str
	  Input sequence

	Returns
	-------
	gcperc : float
	  GC percentage of the input sequence

	Example
	-------
	GetGC("ACAC")
	> 0.5
	"""
	seq = seq.upper()
	STRLen = len(seq)
	CGCount = seq.count('C') + seq.count('G')
	gcperc = CGCount/STRLen		
	return gcperc

def GenerateVariableRegion(chrom, str_start, str_end, \
						alen, ref_motif, motif, maxlen_bp, \
						max_context_size, genome):
	"""
	Main function to generate variable regions
	with STR + context

	Arguments
	---------
	chrom : str
	  Chromosome of the STR
	str_start : int
	  Start coord of the STR
	str_end : int
	  End coord of the STR
	alen : int
	  Allele length to generate (in copy number)
	ref_motif : str
	  Repeat unit sequence in the reference
	motif : str
	  Repeat unit sequence to synthesize
	maxlen_bp : int
	  Maximum allele length we will generate (in bp)
	max_context_size : int
	  Maximum context length (bp)
	genome : pyfaidx.Fasta
	  Reference genome

	Returns
	-------
	var_region : str
	   Variable region
	match_ref : bool
	   True if exactly matches ref
	"""
	match_ref = False
	### Extract flanks and STR seq from reference
	left_flank = str(genome[chrom][str_start-int(np.ceil(max_context_size/2)):str_start]).upper()
	str_refseq = str(genome[chrom][str_start:str_end]).upper()
	right_flank = str(genome[chrom][str_end:str_end+int(np.floor(max_context_size/2))]).upper()

	### Constuct STR region
	str_region = "" # will fill in with cases below

	# Case 1: the motif is the same as the reference. 
	# wnat to keep same imperfections as in the reference
	# to be like the old array
	if motif == "ref":
		ref_alen = int((str_end-str_start+1)/len(ref_motif))
		if alen<ref_alen:
			# deletion, just remove some number of repeat units from the end
			delbp = (ref_alen-alen)*len(ref_motif)
			str_region = str_refseq[0:-1*delbp]
		elif alen == ref_alen:
			# same as reference allele
			str_region = str_refseq
			match_ref = True
		else:
			# insertion: add some number of copies to ref
			ref_motif_rotation = str_refseq[-1*len(ref_motif):] # take from last len(motif) bp of reference seq
			str_region = str_refseq + ref_motif_rotation*(alen-ref_alen)

	# Case 2: the motif is different than the reference
	# Replace with perfect copies
	elif "random" not in motif:
		str_region = motif*alen

	# Case 3: random sequence
	else:
		seq_len = len(ref_motif)*alen
		if motif == "random_matchedGC":
			str_region = GenerateRandomSequence(seq_len, gcperc=GetGC(ref_motif))
		else:
			str_region = GenerateRandomSequence(seq_len, gcperc=0.5)

	### Return entire variable region
	return (left_flank + str_region + right_flank), match_ref

def GetFiller(len_var_reg):
	"""
	Generate filler sequence. Compute length
	based on probelen-(length of all other elements)
	and grab sequence from GLOBAL_FILLER_SEQ

	Arguments
	---------
	len_var_reg : int
	   Length of the variable region

	Returns
	-------
	filler_seq : str
	   Sequence of the filler to use
	"""
	required_elts_len = len(FIVE_PRIME_ADAPT) + len(GIBSON_ASISI) + \
		len(BSAI_RECOG) + len(GIBSON_BSAI_CUT)
	filler_len = PROBE_LEN-len_var_reg-required_elts_len
	filler_seq = GLOBAL_FILLER_SEQ[0:filler_len]
	return filler_seq

def GenerateOligo(vreg, debug=False):
	"""
	Generate the oligo based on the variable region given

	Arguments
	---------
	vreg : str
	   Sequence of the variable region
	debug : bool
	   If true, print out debugging info about the oligo

	Returns
	-------
	oligo_list : list of str
	   List of oligo sequence components to include on the array
	"""
	filler_seq = GetFiller(len(vreg))
	oligo_list = [FIVE_PRIME_ADAPT, vreg, GIBSON_ASISI, filler_seq, BSAI_RECOG, GIBSON_BSAI_CUT]
	assert(len(''.join(oligo_list))==PROBE_LEN)
	return oligo_list

def CheckCutSites(oligo):
	"""
	Check if cutsites to remove the filler sequence are in the oligo
	
	Arguments
	---------
	oligo : str
		Oligo sequence to check for cut sites

	Returns
	-------
	passed : bool
		True if no cut sites found. else False
	"""
	passed = True
	cutsites_to_check = [AsiSI_RECOG, BSAI_RECOG]
	for cutsite in cutsites_to_check:
		if cutsite in oligo:
			passed = False
	return passed

def main():
	### Set up argument parsing ###
	parser = argparse.ArgumentParser(__doc__)
	# Inputs
	parser.add_argument("--reffa", help="Reference genome (fasta)", \
		type=str, required=True)
	parser.add_argument("--rptsbed", help="List of repeats in bed format. "
		"Columns: chrom, start, end, rptunit", \
		type=str, required=True)
	# Outputs
	parser.add_argument("--out", help="Prefix for output files", \
		type=str, required=True)
	# Extra optional arguments
	parser.add_argument("--debug", help="Print helpful debug info", action="store_true")
	## Parse the arguments
	args = parser.parse_args()
	if not os.path.exists(args.reffa):
		sys.stderr.write("Quitting. Can't find " + args.reffa + "!\n")
		sys.exit(1)
	if not os.path.exists(args.rptsbed):
		sys.stderr.write("Quitting. Can't find " + args.rptsbed + "!\n")
		sys.exit(1)
	# Check the arguments
	sys.stderr.write("Building uber array\n")
	sys.stderr.write("Using reference genome: " + args.reffa + "\n")
	sys.stderr.write("Using repeats bed: " + args.rptsbed + "\n")
	sys.stderr.write("Using output prefix: " + args.out + "\n")

	# Set up reference genome. Extract with genome[chr][start:end]
	genome = pyfaidx.Fasta(args.reffa)

	######## First pass ###########
	# Set up output file
	f_oligo = open(args.out + ".oligos.tab", "w")
	f_split = open(args.out + ".oligos.split.tab", "w")

	# Process one STR locus at a time
	with open(args.rptsbed, "r") as f:
		for line in f:
			probecount = 0 # keep track of num probes for this STR
			# Extract relevant info for the STR
			items = line.strip().split()
			chrom = items[0]
			str_start = int(items[1])
			str_end = int(items[2])
			repeat_unit = utils.GetCanonicalMotif(items[3])
			locname="%s_%s_%s_%s"%(chrom, str_start, str_end, repeat_unit)
			if args.debug:
				sys.stderr.write("Processing %s:%s-%s repeat unit=%s\n"%(chrom, str_start, \
						str_end, repeat_unit))

			# Generate the oligos for this STR
			# Planned permutations (initial planning... this is changing slightly as we actually implement it)
			# Different total repeat lengths (homo=10, di=15, tri=25, tetra=20, penta=16, hexa=13)
			# Change motif (6 total, 2 different each of di-, tri-, tetra-)
			# Random sequence of different lengths (2 times each length, matched or varying GC content)
			# Sequence interruptions (will do this later...)
			# Positive and negative strand for each

			# Number of different repeat lengths to generate
			num_rpt_lengths = TOTAL_LENGTHS[len(repeat_unit)]
			allele_lengths = GetAlleleLengths(num_rpt_lengths, len(repeat_unit))
			motifs = ["ref"] + [repeat_unit] + GetAlternateMotifs(exclude=repeat_unit) + ["random"] + ["random_matchedGC"]
			max_motif_len = len([len(m) for m in motifs if "random" not in m])
			max_bp_len = np.min([MAX_STR_LEN, max(allele_lengths)*max_motif_len])

			if utils.ReverseComplement(repeat_unit) != repeat_unit:
				motifs.append(utils.ReverseComplement(repeat_unit))

			### Determine context amount ###
			required_elts_len = len(FIVE_PRIME_ADAPT) + len(GIBSON_ASISI) + \
				len(BSAI_RECOG) + len(GIBSON_BSAI_CUT)
			max_context_size = PROBE_LEN - required_elts_len - max_bp_len
			if args.debug: sys.stderr.write("   [%s] Context size=%s\n"%(locname, max_context_size))

			# Loop through all desired permutations
			for alen in allele_lengths:
				for motif in motifs:
					if args.debug:
						if args.debug: sys.stderr.write("   [%s] Generating oligos with allele len=%s and motif=%s\n"%(locname, alen, motif))
					if "random" not in motif and "ref" not in motif and len(motif)*alen>MAX_STR_LEN:
						if args.debug: sys.stderr.write("   [%s] Skipping len=%s and motif=%s. Too long!\n"%(locname, alen, motif))
						continue
					# Step 1: generate the variable flanking + STR region 
					variable_region, match_ref = GenerateVariableRegion(chrom, str_start, str_end, \
															alen, repeat_unit, motif, \
															max_bp_len, max_context_size, \
															genome)
					# Step 2: Generate the full oligo for both the sequence and the reverse complement
					oligo_num = 0
					for vreg in [variable_region, utils.ReverseComplement(variable_region)]:
						if not CheckCutSites(FIVE_PRIME_ADAPT + vreg):
							if args.debug: sys.stderr.write("   [%s] Skipping len=%s and motif=%s. Failed cut site check.\n"%(locname, alen, motif))
							oligo_num += 1 # still ned to increment this
							continue
						oligo_name = "_".join([chrom, str(str_start), str(str_end), repeat_unit, str(alen), motif, str(oligo_num)])
						if match_ref:
							oligo_name += "*"
						oligo_list = GenerateOligo(vreg, debug=args.debug)
						f_oligo.write("\t".join([oligo_name, ''.join(oligo_list)])+"\n")
						f_split.write("\t".join([oligo_name] + oligo_list)+"\n")
						oligo_num += 1
						probecount += 1
			sys.stderr.write("   [%s] Generated %s probes for %s. Context len=%s\n"%(locname, probecount, locname, max_context_size))
	f_oligo.close()
	f_split.close()

	######## Second pass to remove any redundant probes ###########
	probes = set()
	f_filt = open(args.out + ".oligos.filtered.tab", "w")
	with open(args.out + ".oligos.tab", "r") as f:
		for line in f:
			probe = line.strip().split()[1]
			if probe in probes:
				if args.debug: sys.stderr.write("Removing %s. Redundant probe.\n"%line.strip().split()[0])
			else:
				probes.add(probe)
				f_filt.write(line.strip()+"\n")
	f_filt.close()
	sys.exit(0)

if __name__ == "__main__":
	main()