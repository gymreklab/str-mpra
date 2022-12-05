#!/usr/bin/env python3
"""
Design uber STR array sequence for MPRA

Example command:
./uber-array-design.py \
   --reffa chr21.fa \
   --rptsbed chr21_test.bed \
   --out chr21_test
"""


import argparse
import os
import sys
import pyfaidx
from trtools.utils import utils

# Global variables for oligo construction. These cannot change
FIVE_PRIME_ADAPT = 'ACTGGCCGCTTGACG'
GIBSON_ASISI = 'CACTGCGGCTCCTGCGATCGC'
GLOBAL_FILLER_SEQ = 'TAACCAGGCGTGTTAGCTGCTGTGCTGTCCTACGAGTAAACAGTAGAGTCCGTGGGCACGCGAGCACGGTGAGTCGACTCTGGCCTCATCACCATTTAGTTTGCGCAAGCGCTCTTTTTATAGGACCTGTCTTACATCCCTCATTAACGGAATCGATTACCGGCTAGCGTTGAAATGGAGAAACCGGCTTGCAGTCGAAA'
BSAI_RECOG = 'GGTCTC'
GIBSON_BSAI_CUT = 'TGTCGATCGCGTCGACGAAC'
PROBE_LEN = 230

# Other configuration variables that can be played with
TOTAL_LENGTHS = {1: 10, 2: 15, 3: 25, 4: 20, 5: 16, 6: 13}

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

def GenerateVariableRegion():
	"""
	Main function to generate variable regions
	with STR + context

	TODO: fill out the arguments this will need (ref genome, motif, ...)

	Returns
	-------
	var_region : str
	   Variable region
	"""
	return "XXXXXX" # TODO Not yet implemented

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

def GenerateOligo(vreg):
	"""
	Generate the oligo based on the variable region given

	Arguments
	---------
	vreg : str
	   Sequence of the variable region

	Returns
	-------
	oligo : str
	   Oligo sequence to include on the array
	"""
	filler_seq = GetFiller(len(vreg))
	oligo = FIVE_PRIME_ADAPT + vreg + GIBSON_ASISI + filler_seq + BSAI_RECOG + GIBSON_BSAI_CUT
	return oligo

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

	# Set up output file
	f_oligo = open(args.out + ".oligos.csv", "w")

	# Process one STR locus at a time
	with open(args.rptsbed, "r") as f:
		for line in f:
			# Extract relevant info for the STR
			items = line.strip().split()
			chrom = items[0]
			str_start = int(items[1])
			str_end = int(items[2])
			repeat_unit = utils.GetCanonicalMotif(items[3])

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
			motifs = [repeat_unit] + GetAlternateMotifs(exclude=repeat_unit) + ["random"] + ["random_matchedGC"]
			if utils.ReverseComplement(repeat_unit) != repeat_unit:
				motifs.append(utils.ReverseComplement(repeat_unit))

			# Loop through all desired permutations
			for alen in allele_lengths:
				for motif in motifs:
					if args.debug:
						sys.stderr.write("   Generating oligos with allele len=%s and motif=%s\n"%(alen, motif))
					# Step 1: generate the variable flanking + STR region 
					variable_region = GenerateVariableRegion()
					# Step 2: Generate the full oligo for both the sequence and the reverse complement
					oligo_num = 0
					for vreg in [variable_region, utils.ReverseComplement(variable_region)]:
						oligo_name = "_".join([chrom, str(str_start), str(str_end), repeat_unit, str(alen), motif, str(oligo_num)])
						oligo = GenerateOligo(vreg)
						f_oligo.write(",".join([oligo_name, oligo])+"\n")
						oligo_num += 1

	f_oligo.close()
	sys.exit(0)

if __name__ == "__main__":
	main()