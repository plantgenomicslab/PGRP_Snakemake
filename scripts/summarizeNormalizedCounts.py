#!/usr/bin/env python

# Usage: ./summarizeNormalizedCounts.py [counts_file]
# This script takes a counts file (tsv) and computes gene-wise averages and
# standard deviations among replicates. Treatment/replicate relationships
# are defined by the sraRunsbyExperiment.tsv input file.

import os
import pandas as pd
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sample_table

gene_index = sys.argv[1]
COUNTS_FILE = sys.argv[2]
counts = pd.read_csv(COUNTS_FILE, sep="\t", index_col=gene_index)

SAMPLES_FILE = sample_table.load_sample_table()
# Accepts either 'Treatment' (canonical) or the legacy 'Sample' header — this
# script used to hard-code 'Sample' while the Snakefile hard-coded 'Treatment',
# so one control file could not satisfy both (issue #9).
try:
    GROUP_COLUMN = sample_table.resolve_group_column(SAMPLES_FILE.columns)
except sample_table.SampleTableError as err:
    sys.exit(f"{err}\nExiting...")
REPLICATE_LOOKUP = SAMPLES_FILE.groupby(GROUP_COLUMN)['Replicate'].unique().apply(list).to_dict()

averages = pd.DataFrame()
stdDevs = pd.DataFrame()

for treatment in REPLICATE_LOOKUP:
    replicates = [counts[replicate] for replicate in REPLICATE_LOOKUP[treatment]]
    reps = pd.DataFrame(replicates).transpose()

    averages[treatment] = reps.mean(axis=1)
    stdDevs[treatment] = reps.std(axis=1)

averages.to_csv(COUNTS_FILE + ".average.tsv", sep="\t", index=True)
stdDevs.to_csv(COUNTS_FILE + ".stdDev.tsv", sep="\t", index=True)
