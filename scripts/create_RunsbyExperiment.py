import os
import pandas as pd
import sys
import datetime
import argparse
from itertools import combinations

def create_table(directory):
    df_list = []  # create an empty list to store each row as a DataFrame
    for filename in os.listdir(directory):
        if filename.endswith("R1.fastq.gz") or filename.endswith("R1.fq.gz"):
            run = filename.split('_R1')[0]
            experiment = run.rsplit('_', 1)[0]
            replicate = run.rsplit('_', 1)[1]
            sample = run
            # create a single row DataFrame and append it to the list
#            df_list.append(pd.DataFrame({'Run': [run], 'Experiment': [experiment], 'Replicate': [replicate], 'Sample': [sample]}))
            df_list.append(pd.DataFrame({'Run': [run], 'Experiment': [experiment], 'Replicate': [run], 'Sample': [experiment]}))
    df = pd.concat(df_list, ignore_index=True)  # concatenate all the DataFrames in the list
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script generates RunsByExperiment and **SAMPLE** sample_contrasts files from a given directory.')
    parser.add_argument('directory', help='The directory where the files are located.')
    args = parser.parse_args()
    
    df = create_table(args.directory)
    df.sort_values(by=['Run'], inplace=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    df.to_csv(f'RunsByExperiment_{timestamp}.txt', sep='\t', index=False)

    # Pairwise comparisons for sample contrast
    #pairwise_df = pd.DataFrame(list(combinations(df['experiment'], 2)), columns=['Sample1', 'Sample2'])
    #pairwise_df = pd.DataFrame(list(combinations(df['Experiment'], 2)), header=FALSE)
    #pairwise_df.to_csv(f'output_pairwise_{timestamp}.tsv', sep='\t', index=False)
    pairwise_combinations = list(combinations(df['Experiment'], 2))
    pairwise_combinations = list(set(tuple(sorted(sub)) for sub in pairwise_combinations))
    pairwise_df = pd.DataFrame(pairwise_combinations, columns=['Sample1', 'Sample2'])
    pairwise_df = pairwise_df[pairwise_df['Sample1'] != pairwise_df['Sample2']]
    pairwise_df.sort_values(by=['Sample1', 'Sample2'], inplace=True)
    pairwise_df.to_csv(f'sample_contrasts_{timestamp}.txt', sep='\t', index=False, header=False)
