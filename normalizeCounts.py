#!/usr/bin/env python

import pandas as pd
import sys

mode = sys.argv[1]
GTF_file = sys.argv[2]
input_file = sys.argv[3]
output_prefix = sys.argv[4]

def sumCountLength(gene_id, raw_count, geneLengths):
    sum_per_length = raw_count/geneLengths.loc[gene_id]["Length"]
    return(sum_per_length)

if mode == "HTseq":
    #Load raw counts from HTseq output
    counts = pd.read_csv(input_file, sep="\t", index_col="gene")
    counts = counts[:-5]

    # Compute gene lengths from GTF file
    # These lengths are not congruent with the featureCounts lengths in some cases, need to revisit this
    gtf = pd.read_csv(GTF_file, header=None, sep="\t")
    gtf_exon = gtf.loc[gtf[2] == "exon"]
    gene_length = gtf_exon[8].str.split(";", expand=True)[1]
    gene_length = pd.DataFrame(gene_length)
    gene_length["Length"] = gtf_exon[4] - gtf_exon[3]
    gene_length.columns = ["Geneid", "Length"]
    gene_length["Geneid"] = gene_length["Geneid"].str.replace("[(gene_id )(\")]", "", regex=True)
    gene_length = gene_length.groupby(['Geneid']).sum()

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
