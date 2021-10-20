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
    if "," in counts: return "NOTUNIQUE"
    seq, count = counts.split(":")
    count = int(count)
    if count < 3: return "DEPTH"
    return "PASS"
    
# Filter to require:
# BC/seq assoc seen > 2
# didn't see same BC associated with more than 1 seq
bcsum = data.groupby("barcode", as_index=False).agg({"seq": GetCount})
bcsum["filter"] = bcsum["seq"].apply(GetFilter)

bcsum.to_csv(sys.stdout, sep="\t", index=False)

