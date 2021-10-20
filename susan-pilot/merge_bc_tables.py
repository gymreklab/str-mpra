#!/usr/bin/env python3

import pandas as pd
import sys

fillers = pd.read_csv(sys.argv[1], sep="\t", names=["readname","filler.length","filler"])
barcodes = pd.read_csv(sys.argv[2], sep="\t", names=["readname", "barcode"])
varseq = pd.read_csv(sys.argv[3], sep="\t", names=["readname", "seq"])

fillers["readname"] = fillers["readname"].apply(lambda x: str(x).strip().split()[0][1:])
barcodes["readname"] = barcodes["readname"].apply(lambda x: str(x).strip().split()[0][1:])
varseq["readname"] = varseq["readname"].apply(lambda x: str(x).strip().split()[0])

data = pd.merge(fillers, barcodes, on=["readname"])
data = pd.merge(data, varseq, on=["readname"])

data[["readname","barcode","seq","filler.length"]].to_csv(sys.stdout, index=False, sep="\t")
