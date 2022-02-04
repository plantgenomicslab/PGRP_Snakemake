#!/usr/bin/env python3

# usage: joinSraRelations.py [SRA Project Id]
# This script will download Run and Experiment level metadata from NCBI's SRA
# database and join the information into a single tabele (csv) relating run
# ids (SRR...) to experiment ids (SRX....). 

import pandas as pd
import os
import sys

os.system("curl -o 'sraRunInfo.csv' 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + sys.argv[1] +"'")
os.system("curl -o 'sraExperimentSummary.xml' 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&View=docsumcsv&term="+ sys.argv[1] +"&ContentType=csv&Mode=file'")

runs = pd.read_csv("sraRunInfo.csv")
experiments = pd.read_xml("sraExperimentSummary.xml", xpath="//EXPERIMENT_PACKAGE_SET/EXPERIMENT_PACKAGE/EXPERIMENT")

runs = runs[["Run","Experiment"]]

experiments = experiments[["accession", "TITLE"]]
experiments = experiments.rename(columns={"accession":"Experiment", "TITLE":"Treatment"})
experiments[['alias','Treatment', 'Genus', 'Species', 'Type']] = experiments['Treatment'].astype("string").str.split(' ',expand=True)
experiments = experiments[['Experiment', 'Treatment']]
experiments['Treatment'] = experiments['Treatment'].str.replace(r';', '')

joined = pd.merge(runs, experiments, on="Experiment")

joined.to_csv("sraRunsbyExperiment.csv", index=False)
