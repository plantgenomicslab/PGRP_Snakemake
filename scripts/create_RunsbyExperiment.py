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
            df_list.append(pd.DataFrame({'Run': [run], 'Treatment': [experiment], 'Replicate': [run], 'Sample': [experiment]}))
    if not df_list:
        raise SystemExit(
            f"No R1 read files found in {directory}.\n"
            "This script expects paired-end files named '<sample>_rep<N>_R1.fastq.gz' "
            "or '<sample>_rep<N>_R1.fq.gz'.\n"
            "For other naming schemes use scripts/Experiment_name_composer.py.\n"
            "Exiting..."
        )
    df = pd.concat(df_list, ignore_index=True)  # concatenate all the DataFrames in the list
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script generates RunsByExperiment and **SAMPLE** sample_contrasts files from a given directory.')
    parser.add_argument('directory', help='The directory where the files are located.')
    args = parser.parse_args()
    
    df = create_table(args.directory)
    df.sort_values(by=['Run'], inplace=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    df.to_csv(f'RunsByExperiment_{timestamp}.tsv', sep='\t', index=False)

    # Pairwise comparisons for sample contrast
    #pairwise_df = pd.DataFrame(list(combinations(df['experiment'], 2)), columns=['Sample1', 'Sample2'])
    #pairwise_df = pd.DataFrame(list(combinations(df['Experiment'], 2)), header=FALSE)
    #pairwise_df.to_csv(f'output_pairwise_{timestamp}.tsv', sep='\t', index=False)
    # create_table() emits 'Treatment', never 'Experiment' — pairing off the
    # latter raised KeyError before every run reached this point.
    # combinations() over the *unique* treatments already yields each unordered
    # pair once and never pairs a treatment with itself.
    pairwise_combinations = list(combinations(sorted(df['Treatment'].unique()), 2))
    pairwise_df = pd.DataFrame(pairwise_combinations, columns=['Sample1', 'Sample2'])
    pairwise_df.sort_values(by=['Sample1', 'Sample2'], inplace=True)
    pairwise_df.to_csv(f'sample_contrasts_{timestamp}.tsv', sep='\t', index=False, header=False)
