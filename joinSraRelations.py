#!/usr/bin/env python3

# usage: joinSraRelations.py [SRA Project Id] [replicate regex] [treatment regex]
# example: ./joinSraRelations.py SRP098160 "ZT\d{1,2}_rep\d" "ZT\d{1,2}"
# This script will download Run and Experiment level metadata from NCBI's SRA
# database and join the information into a single tabele (tsv) relating run
# ids (SRR...) to experiment ids (SRX....). For each treatment/replicate a human 
# readable title is parsed out if the run titles using user input reges expressions.

import pandas as pd
import os
import sys

os.system("curl -o 'sraRunInfo.csv' 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + sys.argv[1] +"'")
os.system("curl -o 'sraExperimentSummary.xml' 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&View=docsumcsv&term="+ sys.argv[1] +"&ContentType=csv&Mode=file'")

runs = pd.read_csv("sraRunInfo.csv")
experiments = pd.read_xml("sraExperimentSummary.xml", xpath="//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/EXPERIMENT")

runs = runs[["Run","Experiment"]]

experiments = experiments[["accession", "TITLE"]]
experiments = experiments.rename(columns={"accession":"Experiment", "TITLE":"Replicate"})
experiments['Replicate'] = experiments['Replicate'].str.extract(r'(' + sys.argv[2] + ')')
experiments['Treatment'] = experiments['Replicate'].str.extract(r'(' +sys.argv[3] + ')')

joined = pd.merge(runs, experiments, on="Experiment")

joined.to_csv("sraRunsbyExperiment.tsv", sep="\t", index=False)
