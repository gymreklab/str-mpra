#!/usr/bin/env python3

import re
import sys

infile = sys.argv[1]

def CigarDiff(cigar):
    if "H" in cigar or "S" in cigar: return None
    diff = 0
    for num1, i_or_d, numg2, m in re.findall('(\d+)([ID])(\d+)?([A-Za-z])?', cigar):
        if i_or_d == "I":
            diff += int(num1)
        else: diff -= int(num1)
    return diff

with open(infile, "r") as f:
    for line in f:
        if "UMI" in line: continue
        items = line.strip().split()
        umi, tag, exp, cigar = items
        exp = int(exp)
        obs = CigarDiff(cigar)
        if obs is None:
            diff = None
        else: diff = obs-exp
        output = [umi, tag, exp, cigar, obs, diff]
        sys.stdout.write("\t".join([str(item) for item in output])+"\n")
