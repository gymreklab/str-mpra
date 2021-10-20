#!/usr/bin/env python3

#UMIMAP = "/storage/mlamkin/data/eSTR-data/randomized_barcode_method/results/pSTR1_S1_L001_umi_map.txt"
UMIMAP = "barcodes/barcode_assoc_table.pass.tab"

import sys

def ReverseComplement(seq):
    r"""Get reverse complement of a sequence.
    Converts everything to uppsercase.
    Parameters
    ----------
    seq : str
          String of nucleotides.
    Returns
    -------
    revcompseq : str
          Reverse complement of the input sequence.
    Examples
    --------
    >>> ReverseComplement("AGGCT")
    "AGCCT"
    """
    seq = seq.upper()
    newseq = ""
    size = len(seq)
    for i in range(len(seq)):
        char = seq[len(seq)-i-1]
        if char == "A":
            newseq += "T"
        elif char == "G":
            newseq += "C"
        elif char == "C":
            newseq += "G"
        elif char == "T":
            newseq += "A"
        else: newseq += "N"
    return newseq

filelist = sys.argv[1:]

# Get all the counts in a big dictionary
count_dict = {} # UMI -> prefix -> count
labels = []
for fn in filelist:
    label = fn.strip(".counts")
    labels.append(label)
    with open(fn, "r") as f:
        for line in f:
            umi, count = line.strip().split()
            count = int(count)
            if umi not in count_dict: count_dict[umi] = {}
            count_dict[umi][label] = count

# Load map of UMI -> seqname
umi2seq = {}
with open(UMIMAP, "r") as f:
    for line in f:
        umi, seqname = line.strip().split()
        umi2seq[umi] = seqname

sys.stdout.write("\t".join(["UMI", "seqname"]+labels)+"\n")
for umi in count_dict.keys():
#    items = [umi, umi2seq.get(ReverseComplement(umi), "NA")]
    items = [umi, umi2seq.get(umi, "NA")]
    for label in labels:
        items.append(count_dict[umi].get(label, 0))
    sys.stdout.write("\t".join([str(item) for item in items])+"\n")
