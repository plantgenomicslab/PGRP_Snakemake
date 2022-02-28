#!/usr/bin/env python

# Usage: ./mergeTPMCalculator.py [TPMCalculator .out files]
# This script takes multiple TPMCalculator .out files and returns
# a tsv file of TPM values for each file.

import pandas as pd
import sys

files = sys.argv[2:]

output = pd.read_csv(sys.argv[1], sep="\t", index_col=False)[["Gene_Id", "ExonLength", "ExonTPM"]]
output = output.rename(columns={"ExonTPM":sys.argv[1].replace("_genes.out", "").replace("output/counts/tpmcalculator/", "")})

for tpmCalcOut in files:
    cur_out = pd.read_csv(tpmCalcOut, sep="\t", index_col=False)[["Gene_Id", "ExonTPM"]]
    print(cur_out.head())
    cur_out = cur_out.rename(columns={"ExonTPM": tpmCalcOut.replace("_genes.out", "").replace("output/counts/tpmcalculator/", "")})
    output = pd.merge(output, cur_out, on="Gene_Id")

output.to_csv("output/counts/tpmcalculator/tpmcalculator-merged.tsv", sep="\t", index=False)
