import os
import pandas as pd
import sys
import datetime
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sample_table

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
    parser = argparse.ArgumentParser(
        description='This script generates RunsByExperiment and **SAMPLE** sample_contrasts files from a given directory.',
        epilog='Each sample_contrasts line "A<TAB>B" yields log2FoldChange = log2(A/B) in the '
               'Trinity/DESeq2 output, so the treatment belongs first and the control second. '
               'Without --reference pairs are written in alphabetical order; check the '
               'direction and edit the file before running DEG analysis.')
    parser.add_argument('directory', help='The directory where the files are located.')
    parser.add_argument('--reference', metavar='TREATMENT',
                        help='Control/reference treatment: written second in every pair '
                             'that contains it, so its fold changes read "vs TREATMENT".')
    args = parser.parse_args()
    
    df = create_table(args.directory)
    df.sort_values(by=['Run'], inplace=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    df.to_csv(f'RunsByExperiment_{timestamp}.tsv', sep='\t', index=False)

    # Pairwise comparisons for sample contrast.
    # Trinity's DESeq2 contrast for a line 'A<TAB>B' is log2(A/B); see
    # sample_table.pairwise_contrasts for how the direction is chosen.
    try:
        pairwise_combinations = sample_table.pairwise_contrasts(df['Treatment'], args.reference)
    except sample_table.SampleTableError as exc:
        sys.exit(f"ERROR: {exc}")
    pairwise_df = pd.DataFrame(pairwise_combinations, columns=['Sample1', 'Sample2'])
    contrasts_path = f'sample_contrasts_{timestamp}.tsv'
    pairwise_df.to_csv(contrasts_path, sep='\t', index=False, header=False)
    if args.reference is None:
        print(f"NOTE: {contrasts_path} lists pairs alphabetically. Line 'A<TAB>B' gives "
              "log2FC = log2(A/B): put the treatment first, or re-run with --reference <control>.",
              file=sys.stderr)
