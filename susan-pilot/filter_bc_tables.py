#!/usr/bin/env python3

import pandas as pd
import sys

data = pd.read_csv(sys.argv[1], sep="\t")

def GetCount(seqs):
    counts = {}
    for s in seqs:
        counts[s] = counts.get(s, 0) + 1
    return ",".join(["%s:%s"%(s, counts[s]) for s in counts.keys()])

def GetFilter(counts):
    total = 0
    maxcount = 0
    best = None
    #if "," in counts: return "NOTUNIQUE"
    for item in counts.split(","):
        seq, count = item.split(":")
        count = int(count)
        total += count
        if count > maxcount:
            maxcount = count
            best = seq
    if maxcount < 2: return "DEPTH"
    if maxcount*1.0/total < 0.75: return "NOTUNIQUE"
    return "PASS"

def GetBest(counts):
    total = 0
    maxcount = 0
    best = None
    #if "," in counts: return "NOTUNIQUE"
    for item in counts.split(","):
        seq, count = item.split(":")
        count = int(count)
        total += count
        if count > maxcount:
            maxcount = count
            best = seq
    return best

# Require filler sequence matches expected
def GetExpFiller(seq):
    if seq == "*": return -1
    if seq == "1_POMC_-22": return 104
    if seq == "2_POMC_0": return 82
    if seq == "3_POMC_9": return 73
    if seq == "4_POMC_20": return 64
    if seq == "5_POMC_29": return 55
    return 0
data["exp.filler"] = data.apply(lambda x: GetExpFiller(x["seq"]), 1)

sys.stderr.write("Before filter by filler seq match var seq: %s\n"%data.shape[0])
data = data[data["exp.filler"] == data["filler.length"]]
sys.stderr.write("After filter by filler seq match var seq: %s\n"%data.shape[0])
    
# Filter to require:
# BC/seq assoc seen >= 2
# didn't see same BC associated with more than 1 seq
bcsum = data.groupby("barcode", as_index=False).agg({"seq": GetCount})
bcsum["filter"] = bcsum["seq"].apply(GetFilter)
bcsum["orig.seq"] = bcsum["seq"]
bcsum["seq"] = bcsum["seq"].apply(GetBest)

bcsum.to_csv(sys.stdout, sep="\t", index=False)

