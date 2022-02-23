#!/usr/bin/env python

# Usage: ./summarizeNormalizedCounts.py [counts_file]
# This script takes a counts file (tsv) and computes gene-wise means and
# standard deviations among replicates. Treatment/replicate relationships
# are defined by the sraRunsbyExperiment.tsv input file.

import pandas as pd
import sys

COUNTS_FILE = sys.argv[1]
counts = pd.read_csv(COUNTS_FILE, sep="\t", index_col="Geneid")

SAMPLES_FILE = pd.read_csv("sraRunsbyExperiment.tsv", sep="\t")
REPLICATE_LOOKUP = SAMPLES_FILE.groupby("Treatment")['Replicate'].unique().apply(list).to_dict()

averages = pd.DataFrame()
stdDevs = pd.DataFrame()

for treatment in REPLICATE_LOOKUP:
    replicates = [counts[replicate] for replicate in REPLICATE_LOOKUP[treatment]]
    reps = pd.DataFrame(replicates).transpose()

    averages[treatment] = reps.mean(axis=1)
    stdDevs[treatment] = reps.std(axis=1)

averages.to_csv(COUNTS_FILE + ".mean.tsv", sep="\t", index=True)
stdDevs.to_csv(COUNTS_FILE + ".stdDev.tsv", sep="\t", index=True)
