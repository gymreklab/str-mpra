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

	# Process one STR locus at a time
	with open(args.rptsbed, "r") as f:
		for line in f:
			# Extract relevant info for the STR
			items = line.strip().split()
			chrom = items[0]
			str_start = int(items[1])
			str_end = int(items[2])
			repeat_unit = items[3]

			if args.debug:
				sys.stderr.write("Processing %s:%s-%s repeat unit=%s\n"%(chrom, str_start, \
						str_end, repeat_unit))

	sys.exit(0)

main()