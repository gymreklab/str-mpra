#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

data = pd.read_csv(sys.argv[1], sep="\t")

# Filter if no seqname assigned
data = data[~(data["seqname"].apply(str) == "nan")]

# SB_pSTR1_DNA1_S19  SB_pSTR1_DNA2_S20  SB_pSTR1_DNAL_S21  SB_pSTR1_RNA1_S16  SB_pSTR1_RNA2_S17  SB_pSTR1_RNAL_S18
# Get RNA/DNA ratio for each condition
def GetRatio(dnacol, rnacol):
    if dnacol < 10 or rnacol < 10: return np.nan
    else: return rnacol*1.0/dnacol

for transfections in [("SB_pSTR1_DNA1_S19", "SB_pSTR1_RNA1_S16", "1"), ("SB_pSTR1_DNA2_S20", "SB_pSTR1_RNA2_S17", "2"), ("SB_pSTR1_DNAL_S21", "SB_pSTR1_RNAL_S18", "L")]:
    dnacol, rnacol, label = transfections
    data["ratio.%s"%label] = data.apply(lambda x: GetRatio(x[dnacol], x[rnacol]), 1)

# Filter if all are nan
data = data[~(np.isnan(data["ratio.1"]) & (np.isnan(data["ratio.2"])) & (np.isnan(data["ratio.L"])))]
data = data.sort_values("seqname")

# Print counts to a file
data.to_csv(sys.argv[2], sep="\t", index=False)

# Get a summary
print(data.groupby("seqname", as_index=False).agg({"ratio.1": np.mean, "ratio.2": np.mean, "ratio.L": np.mean, "UMI": len}))

