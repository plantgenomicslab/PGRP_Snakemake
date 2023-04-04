#!/usr/bin/env python
##  ls -1 /data/gpfs/assoc/pgl/data/iceplant/geneatlas_05042022/raw/*CL* | sed 's/.*\///g; s/\.R.*//g' | sort -u > CL_list
## 
## python script.py input.txt RunsByExperiment.tsv

import sys
import os

# get command line arguments
script, input_file, output_file = sys.argv

# define input data
with open(input_file, 'r') as f:
    input_data = f.read().splitlines()

# define output file headers
headers = ['Run', 'Replicate', 'Sample']

# define function to extract sample name
def extract_sample_name(input_str):
    # remove _rep* or _Rep* from the end of the input string
    sample_name = input_str.rsplit('_', 1)[0]
    return sample_name

# define function to create output rows
def create_output_rows(input_data):
    output_rows = []
    for input_str in input_data:
        sample_name = extract_sample_name(input_str)
        output_rows.append([input_str, input_str, sample_name])
    return output_rows

# create output rows
output_rows = create_output_rows(input_data)

# write output file
with open(output_file, 'w') as f:
    # write headers
    f.write('\t'.join(headers) + '\n')
    # write rows
    for row in output_rows:
        f.write('\t'.join(row) + '\n')

