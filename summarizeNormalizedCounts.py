#!/usr/bin/env python

# Usage: ./summarizeNormalizedCounts.py [counts_file]
# This script takes a counts file (tsv) and computes gene-wise averages and
# standard deviations among replicates. Treatment/replicate relationships
# are defined by the sraRunsbyExperiment.tsv input file.

import pandas as pd
import sys

gene_index = sys.argv[1]
COUNTS_FILE = sys.argv[2]
counts = pd.read_csv(COUNTS_FILE, sep="\t", index_col=gene_index)

SAMPLES_FILE = pd.read_csv("sraRunsbyExperiment.tsv", sep="\t")
REPLICATE_LOOKUP = SAMPLES_FILE.groupby("Treatment")['Replicate'].unique().apply(list).to_dict()

averages = pd.DataFrame()
stdDevs = pd.DataFrame()

for treatment in REPLICATE_LOOKUP:
    replicates = [counts[replicate] for replicate in REPLICATE_LOOKUP[treatment]]
    reps = pd.DataFrame(replicates).transpose()

    averages[treatment] = reps.mean(axis=1)
    stdDevs[treatment] = reps.std(axis=1)

averages.to_csv(COUNTS_FILE + ".average.tsv", sep="\t", index=True)
stdDevs.to_csv(COUNTS_FILE + ".stdDev.tsv", sep="\t", index=True)
