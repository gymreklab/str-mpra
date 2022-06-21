#!/usr/bin/env python3
"""
./annoate_bam.py inbam outbam
"""

import pysam
import sys
from Bio.Seq import Seq

TAGINFOFILE = "/storage/mlamkin/projects/eSTR-data/labels_for_bam_analysis.txt"

try:
    inbam_f = sys.argv[1]
    outbam_f = sys.argv[2]
except:
    sys.stderr.write(__doc__)

# Load tag info (tag -> count). Include forward and reverse complement
def ReverseComplement(seq):
    fseq = Seq(seq)
    return str(fseq.reverse_complement())

taginfo = {}
with open(TAGINFOFILE, "r") as f:
    for line in f:
        if "Start" in line: continue
        if "Fun" in line: continue
        items = line.strip().split()
        tag = items[7]
        num = items[6]
#        taginfo[tag] = num
        taginfo[ReverseComplement(tag)] = num

# Set up bam readers/writers
inbam = pysam.AlignmentFile(inbam_f, "rb")
outbam = pysam.AlignmentFile(outbam_f, "wb", template=inbam)

MAXCLIP = 5
# Iterate through reads
for read in inbam.fetch():
    # Filter if lots of clipping
    filt = False
    for item in read.cigartuples:
        if (item[0] in [5,4]) and (item[1] >= MAXCLIP):
            filt=True
            break
    if filt: continue
    # Annotate tag
    read_tag = read.query_name.split("_")[-1]
    count = taginfo.get(read_tag, None)
    if count is None: continue
    read.tags += [("RN",count)]
    outbam.write(read)

inbam.close()
outbam.close()
