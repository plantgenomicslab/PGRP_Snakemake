#!/usr/bin/env python3

# Usage: ./makeRSEMMatrix.py RunsByExperiment.tsv inputDir column

import pandas as pd
import sys

if len(sys.argv) < 4:
    sys.exit("Not enough arguments...\nUsage: ./makeRSEMMatrix.py RunsByExperiment.tsv inputDir column")

#load RunsByExperiment.tsv and extract replicate names
rbe = pd.read_csv(sys.argv[1], sep="\t")
rbe = list(rbe["Replicates"].values)

matrix = pd.DataFrame()

# Read in RSEM results for each replicate and join on gene ids
for i in range(len(rbe)):
    RSEM_result = pd.read_csv(f"{sys.argv[2]}/{rbe[i]}.genes.results", sep="\t")
    
    if i == 0:
        matrix["gene_id"] = RSEM_result["gene_id"]
    
    RSEM_result = RSEM_result[["gene_id", sys.argv[3]]]
    RSEM_result.columns = ["gene_id", rbe[i]]
    matrix=matrix.merge(RSEM_result, on="gene_id", how="left")

matrix.to_csv(f"{sys.argv[2]}/RSEM_{sys.argv[3]}.tsv", index=None, sep="\t")

