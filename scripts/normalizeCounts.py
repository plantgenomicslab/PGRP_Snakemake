#!/usr/bin/env python

import pandas as pd
import sys

mode = sys.argv[1]
GTF_file = sys.argv[2]
input_file = sys.argv[3]
output_prefix = sys.argv[4]
ref_dir = sys.argv[5]

def sumCountLength(gene_id, raw_count, geneLengths):
    sum_per_length = raw_count/geneLengths.loc[gene_id]["Length"]
    return(sum_per_length)

if mode == "HTseq":
    #Load raw counts from HTseq output
    counts = pd.read_csv(input_file, sep="\t", index_col="gene")
    counts = counts[:-5]
    
    # Load exon lengths
    gene_length = pd.read_csv(ref_dir + "cds_length.tsv", header=None, sep="\t")
    gene_length.columns = ["Geneid", "Length"]
    gene_length = gene_length.set_index("Geneid")
    print(gene_length.head())

    #Compute FPKM
    sum_count = counts.sum()
    fpkm = counts.apply(lambda row: (row/(gene_length.loc[row.name]["Length"]*sum_count))*10**9, axis=1)

    #Compute TPM
    sum_count_length = counts.apply(lambda row: sumCountLength(row.name, row, gene_length) , axis=1)
    sum_count_length = sum_count_length.sum()
    tpm = counts.apply(lambda row : (row/(gene_length.loc[row.name]["Length"]*sum_count_length))*10**6 , axis=1)

elif mode == "featureCounts":
    #Load raw counts from featureCounts output
    counts = pd.read_csv(input_file, sep="\t", index_col="Geneid")
    counts_only = counts.iloc[:, 5:]

    #Compute FPKM
    sum_count = counts_only.sum()
    fpkm = counts_only.apply(lambda row: (row/(counts.loc[row.name]["Length"]*sum_count))*10**9, axis=1)

    #Compute TPM
    sum_count_length = counts_only.apply(lambda row: sumCountLength(row.name, row, counts) , axis=1)
    sum_count_length = sum_count_length.sum()
    tpm = counts_only.apply(lambda row: (row/(counts.loc[row.name]["Length"]*sum_count_length))*10**6, axis=1)

tpm.to_csv(output_prefix + ".tpm.tsv", sep="\t", index=True)
fpkm.to_csv(output_prefix + ".fpkm.tsv", sep="\t", index=True)
