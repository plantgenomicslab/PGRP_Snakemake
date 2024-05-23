import sys
import os
import glob
import re
def show_help():
	print("Usage: python Experiment_composer.py <input_folder> [output_file]")
	print("Generates 'RunsByExperiment.tsv' a tab-separated file with columns: Run, Replicate, Sample")
	print("If no output_file is provided, the default 'RunsByExperiment.tsv' will be used.")


# get command line arguments
if len(sys.argv) > 2:
	input_folder = sys.argv[1]
	output_file = sys.argv[2]
else:
	input_folder = sys.argv[1]
	output_file = "RunsByExperiment.tsv"

# define function to generate input data
def generate_input_data(input_folder):
	input_files = []
	for ext in ['.fq', '.fastq', '.fastq.gz', '.fq.gz']:
		input_files.extend(glob.glob(os.path.join(input_folder, f'*{ext}')))

	input_data = []
	for input_file in input_files:
		base_name = os.path.basename(input_file)
		sample_name = base_name.split(".")[0]
		sample_name = sample_name.rsplit("_", 1)[0]
		input_data.append(sample_name)

	input_data = sorted(set(input_data))
	return input_data

# generate input data
input_data = generate_input_data(input_folder)

# define output file headers
headers = ['Run', 'Replicate', 'Treatment']

# define function to create output rows
def create_output_rows(input_data):
	output_rows = []
	rep_pattern = re.compile(r'_rep\d+$|_Rep\d+$')
	for input_str in input_data:
		sample_name = rep_pattern.sub('', input_str)
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
