import yaml, sys, os
import pandas as pd

# Load run and sample information from the sraRunsbyExperiment.tsv input file (user must provide)
SAMPLES_FILE = pd.read_csv("RunsbyExperiment.tsv", sep="\t")
SAMPLE_LIST = list(set(SAMPLES_FILE["Run"].values.tolist()))
REPLICATE_LIST = list(set(SAMPLES_FILE["Replicate"].values.tolist()))
REPLICATE_LOOKUP = SAMPLES_FILE.groupby("Replicate")['Run'].apply(list).to_dict()

# Load configs
configfile: "config.yml"
try:
	with open('config.yaml', "r") as config_file:
		config_dict = yaml.safe_load(config_file)
except yaml.YAMLError as e:
	logging.error("Could not load config file! Check config.yml ... see error below")
	sys.exit(e)

LAYOUT = config["libraryType"] #TODO enable mixed library types
if LAYOUT == "PAIRED":
	PAIR_LIST = config["PAIR_LIST"]

# Check for an indexed reference genome or prepared genome for RSEM
if "RSEM" in config_dict["readCounting"]:
	if not os.path.exists(f"{config_dict['RSEM_prepared_genome']}.seq"):
		sys.exit(f"Cannot locate '{config_dict['RSEM_prepared_genome']}'.\nProvide a reference genome prepared with 'rsem-prepare-reference' to use RSEM.\nExiting...")

for ref in config_dict["ref"]:
	if not os.path.exists(f"{config_dict['genomeDir']}/{ref}"):
		sys.exit("The " + ref + " file is missing.\nHave you provided the correct paths to a reference genome indexed by STAR?\nExiting...")

# Configure DEG calculation (user must provide)
#TODO: write better method for identifying DEG contrasts output files. Maybe 'checkpoint'?
contrasts = os.path.exists(config_dict["sample_contrast"])
if config_dict["runDEG"] == "yes":
	try:
		CONTRASTS_FILE = pd.read_csv(config_dict["sample_contrast"], sep="\t", header=None)
		cTop = CONTRASTS_FILE.iloc[-1][0]
		cBottom = CONTRASTS_FILE.iloc[-1][1]
	except (FileNotFoundError, pd.errors.EmptyDataError):
		cTop = ""
		cBottom = ""
else:
	cTop = ""
	cBottom = ""

# Create sample output folders
os.makedirs("output/logs/", exist_ok=True)
os.makedirs("output/DEG/", exist_ok=True)

counting_options  = ("RSEM","featureCounts","HTseq","TPMcalculator")
for opt in counting_options:
	if opt in config_dict["readCounting"]:
		os.makedirs(f"output/counts/{opt}/", exist_ok=True)

for rep in REPLICATE_LIST:
	os.makedirs("output/" + rep + "/bam/", exist_ok=True)
	for r in REPLICATE_LOOKUP[rep]:
		os.makedirs("output/" + rep + "/" + r  + "/raw/", exist_ok=True)
		os.makedirs("output/" + rep + "/" + r  + "/bam/", exist_ok=True)
		os.makedirs("output/" + rep + "/" + r  + "/logs/", exist_ok=True)
		os.makedirs("output/" + rep + "/" + r  + "/trim/", exist_ok=True)

# Generate a list of all target output files
def allInput():
	inputs = []

	# Count and DEG files
	if "RSEM" in config_dict["readCounting"]:
		inputs += ["output/counts/RSEM/RSEM_TPM.tsv.average.tsv",
					"output/counts/RSEM/RSEM_TPM.tsv"]
		if config_dict["runDEG"] == "yes":
			inputs += ["output/DEG/RSEM_expected_count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results",
						"output/DEG/RSEM_expected_count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset"]
	if "featureCounts" in config_dict["readCounting"]:
		inputs += ["output/counts/featureCounts/featureCounts.cnt",
					"output/counts/featureCounts/featureCounts.tpm.tsv",
					"output/counts/featureCounts/featureCounts.fpkm.tsv",
					"output/counts/featureCounts/featureCounts.tpm.tsv.average.tsv",
					"output/counts/featureCounts/featureCount_clean.cnt"]
		if config_dict["runDEG"] == "yes":
			inputs += ["output/DEG/featureCount_clean.cnt." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset",
						"output/DEG/featureCount_clean.cnt." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results"]
	if "HTseq" in config_dict["readCounting"]:
		inputs += ["output/counts/htseq/htseq-count.tsv",
					"output/counts/htseq/htseq-count.tpm.tsv",
					"output/counts/htseq/htseq-count.fpkm.tsv",
					"output/counts/htseq/htseq-count.tpm.tsv.average.tsv"]
		if config_dict["runDEG"] == "yes":
			inputs += ["output/DEG/htseq-count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results",
						"output/DEG/htseq-count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset"]
	if "TPMcalculator" in config_dict["readCounting"]:
		inputs += ["output/counts/tpmcalculator/tpmcalculator-merged.tsv",
					"output/counts/tpmcalculator/tpmcalculator-merged.tsv.average.tsv"]
		
	# Per-replicate data and processing files
	for replicate in REPLICATE_LIST:
		inputs.append("output/" + replicate + "/bam/" + replicate + ".bam")
		
		if "TPMcalculator" in config_dict["readCounting"]:
			inputs.append("output/counts/tpmcalculator/" + replicate + "_genes.out")
		
		for sample in REPLICATE_LOOKUP[replicate]:
			inputs.append("output/" + replicate + "/" + sample + "/bam/" + sample + ".bamAligned.sortedByCoord.out.bam")
			if LAYOUT == "PAIRED":
				for pair in PAIR_LIST:
					inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample  + pair + ".fastq.gz")
					inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample  + pair + "_fastqc.zip")
					inputs.append("output/" + replicate + "/" + sample + "/trim/" + sample + pair + ".fq.gz")
			elif LAYOUT == "SINGLE":
				inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample + ".fastq.gz")
				inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample + "_fastqc.zip")
				inputs.append("output/" + replicate + "/" + sample + "/trim/" + sample + "_trimmed.fq.gz")
			else:
				print("ERROR: libraryType must be 'SINGLE' or 'PAIRED' ... see config/config.json to set libraryType")
				quit()
	return(inputs)

if LAYOUT == "PAIRED":
	ruleorder: align_PAIRED > align_SINGLE
elif LAYOUT == "SINGLE":
	ruleorder: align_SINGLE > align_PAIRED

localrules: all, importRaw_PAIRED, importRaw_SINGLE

rule all:
	input:
		allInput()

rule importRaw_PAIRED:
	output: "output/{replicate}/raw/{replicate}{pair}.fastq.gz"
	message: "Importing raw data: ln " + config["rawInputDir"] + "/{wildcards.replicate}_{wildcards.pair}.fastq.gz {output}"
	run:
		shell("ln -s " + config["rawInputDir"] + "/{wildcards.replicate}{wildcards.pair}.fastq.gz {output}")

rule importRaw_SINGLE:
	output: "output/{replicate}/raw/{replicate}.fastq.gz"
	message: "Importing raw data: ln " + config["rawInputDir"] + "/{wildcards.replicate}.fastq.gz {output}"
	run:
		shell("ln -s " + config["rawInputDir"] + "/{wildcards.replicate}.fastq.gz {output}")

rule fetchSRA:
	output: config["sra_dir"] + "/{sample}.sra" #"output/sra/{sample}.sra"
	message: "-----Fetching {wildcards.sample} SRA files-----"
	threads: config["threads"]["fetchSRA"]
	run:
		# Download SRA files from NCBI SRA
		shell("prefetch {wildcards.sample} --output-file {output}")
		# Check sra files to make sure they are valid
		shell("echo '--------Validating {wildcards.sample}.sra--------'")
		shell("set +e")
		shell("if ! vdb-validate {output}; then echo 'vdb-validate found errors in {output}. Check log here:{log}. Exiting......' && exit 1; fi")

rule convertSRAtoFastq_PAIRED:
	input: config["sra_dir"] + "/{sample}.sra"
	output: 
		"output/{replicate}/{sample}/raw/{sample}{pair}.fastq.gz",
		"output/{replicate}/{sample}/raw/{sample}{pair}.fastq.gz"
	message: "-----Converting {wildcards.sample} SRA to Fastq files-----"
	threads: config["threads"]["convertSRAtoFastq"]
	log: "output/{replicate}/{sample}/logs/{sample}_fastqdump.log"
	run:
		# Convert sra files to fastQ format
		shell("parallel-fastq-dump --sra-id {wildcards.sample} \
				--threads {threads} --split-e --gzip \
				--outdir output/{wildcards.replicate}/{wildcards.sample}/raw \
				2> {log}")

rule convertSRAtoFastq_SINGLE:
	input: config["sra_dir"] +"/{sample}.sra"
	output: "output/{replicate}/{sample}/raw/{sample}.fastq.gz"
	message: "-----Converting {wildcards.sample} SRA to Fastq files-----"
	threads: config["threads"]["convertSRAtoFastq"]
	log: "output/{replicate}/{sample}/logs/{sample}_fastqdump.log"
	run:
		# Convert sra files to fastQ format
		shell("parallel-fastq-dump --sra-id {wildcards.sample} \
			   --threads {threads} --split-e --gzip \
			   --outdir output/{wildcards.replicate}/{wildcards.sample}/raw \
			   2> {log}")

rule fastqc_raw_PAIRED:
	input:
		"output/{replicate}/{sample}/raw/{sample}{pair}.fastq.gz",
		"output/{replicate}/{sample}/raw/{sample}{pair}.fastq.gz"
	output:
		"output/{replicate}/{sample}/raw/{sample}{pair}_fastqc.zip",
		"output/{replicate}/{sample}/raw/{sample}{pair}_fastqc.zip"
	message: "-----Running Fastqc_raw {wildcards.sample}-----"
	threads: config["threads"]["fastqc_raw"]
	log: "output/{replicate}/{sample}/logs/{sample}_raw_fastqc.log"
	run:
		# Run quality control on raw reads
		shell("fastqc --threads {threads} {input} 2> {log}")

rule fastqc_raw_SINGLE:
	input: "output/{replicate}/{sample}/raw/{sample}.fastq.gz"
	output: "output/{replicate}/{sample}/raw/{sample}_fastqc.zip"
	message: "-----Running Fastqc_raw {wildcards.sample}-----"
	threads: config["threads"]["fastqc_raw"]
	log: "output/{replicate}/{sample}/logs/{sample}_raw_fastqc.log"
	run:
		# Run quality control on raw reads
		shell("fastqc --threads {threads} {input} 2> {log}")

rule trim_PAIRED:
	input:
		fwd_fastq = "output/{replicate}/{sample}/raw/{sample}{pair}.fastq.gz",
		rev_fastq = "output/{replicate}/{sample}/raw/{sample}{pair}.fastq.gz"
	output:
		"output/{replicate}/{sample}/trim/{sample}{pair}.fq.gz",
		"output/{replicate}/{sample}/trim/{sample}{pair}.fq.gz"
	message: "-----Trimming {wildcards.sample}-----"
	log: "output/{replicate}/{sample}/logs/{sample}_trim.log"
	threads: config["threads"]["trim"]
	run:
		# Trim reads using TrimGalore
		shell("trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
			 --cores {threads} --max_n 40 --gzip -o output/{wildcards.replicate}/{wildcards.sample}/trim {input.fwd_fastq} {input.rev_fastq} \
			 2> {log}")
		# Update filenames to improve MultiQC output format
		shell("for file in output/{wildcards.replicate}/{wildcards.sample}/trim/*_val_*; do mv $file $(echo $file | sed s/_val_[0-9]//); done")
		# Run quality control on trimmed reads
		shell("fastqc --threads {threads} {output}")

rule trim_SINGLE:
	input: "output/{replicate}/{sample}/raw/{sample}.fastq.gz"
	output:"output/{replicate}/{sample}/trim/{sample}_trimmed.fq.gz"
	message: "-----Trimming {wildcards.sample}-----"
	log: "output/{replicate}/{sample}/logs/{sample}_trim.log"
	threads: config["threads"]["trim"]
	run:
		# Trim reads using TrimGalore
		shell("trim_galore --three_prime_clip_R1 5 \
			   --cores {threads} --max_n 40 --gzip -o output/{wildcards.replicate}/{wildcards.sample}/trim {input} \
			   2> {log}")
		# Run quality control on trimmed reads
		shell("fastqc --threads {threads} {output}")

rule align_PAIRED:
	input:
		fwd_fastq = "output/{replicate}/{sample}/trim/{sample}{pair}.fq.gz",
		rev_fastq = "output/{replicate}/{sample}/trim/{sample}{pair}.fq.gz"
	output:"output/{replicate}/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam"
	message: "-----Aligning {wildcards.sample}-----"
	log: "output/{replicate}/{sample}/logs/{sample}_align.log"
	threads: config["threads"]["align"]
	run:
		# Calculate run-level alignments to reference
		shell("STAR --runMode alignReads --runThreadN {threads} \
				--outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--genomeDir " + config["genomeDir"]  + " \
			--readFilesCommand gunzip -c --readFilesIn {input.fwd_fastq} {input.rev_fastq} \
				--outSAMtype BAM SortedByCoordinate --outFileNamePrefix output/{wildcards.replicate}/{wildcards.sample}/bam/{wildcards.sample}.bam  \
		2> {log}")
		# Perform check on output bam file to ensure it is not corrupted
		shell("echo '--------Checking {output}----------'")
		shell("set +e")
		shell("if ! samtools flagstat {output}; then echo 'samtools flagstat found errors in {output}. Check log here: {log}. Exiting......' && exit 1; fi")

rule align_SINGLE:
	input: "output/{replicate}/{sample}/trim/{sample}_trimmed.fq.gz"
	output:"output/{replicate}/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam"
	message: "-----Aligning {wildcards.sample}-----"
	log: "output/{replicate}/{sample}/logs/{sample}_align.log"
	threads: config["threads"]["align"]
	run:
		# Calculate run-level alignments to reference
		shell("STAR --runMode alignReads --runThreadN {threads} \
			   --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 \
		   --quantMode TranscriptomeSAM GeneCounts \
			   --genomeDir " + config["genomeDir"]  + " \
			   --readFilesCommand gunzip -c --readFilesIn {input} \
			   --outSAMtype BAM SortedByCoordinate --outFileNamePrefix output/{wildcards.replicate}/{wildcards.sample}/bam/{wildcards.sample}.bam  \
			   2> {log}")
		# Perform check on output bam file to ensure it is not corrupted
		shell("echo '--------Checking {output}----------'")
		shell("set +e")
		shell("if ! samtools flagstat {output}; then echo 'samtools flagstat found errors in {output}. Check log here: {log}. Exiting......' && exit 1; fi")

def mergeInput(wildcards):
	inputs = expand("output/" + wildcards.replicate + "/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=REPLICATE_LOOKUP[wildcards.replicate])
	return(inputs)

rule merge:
	input: mergeInput
	output: "output/{replicate}/bam/{replicate}.bam"
	message: "-----Merging bam files by experiment-----"
	log: "output/{replicate}/bam/{replicate}_merge.log"
	threads: config["threads"]["merge"]
	run:
		# note: bam files must be sorted in same manner prior to merge
		if len(REPLICATE_LOOKUP[wildcards.replicate]) > 1:
			# Merge bam files from same sample
			shell("samtools merge --threads {threads} -o {output} {input}") 
			# Create bam index folder for HTseq
			shell("samtools index {output}") 
		else:
			# Only 1 run in a sample, add it to sample-level output
			shell("cp {input} {output}")
			print("1 run per sample...skipping BAM merge")
		# Perform check on output bam file to ensure it is not corrupted
		shell("echo '--------Checking {output}----------'")
		shell("set +e")
		shell("if ! samtools flagstat {output}; then echo 'samtools flagstat found errors in {output}. Check log here: {log}. Exiting......' && exit 1; fi")

rule TPMCalculator:
	input: "output/{replicate}/bam/{replicate}.bam"
	output: "output/counts/tpmcalculator/{replicate}_genes.out"
	message: "------- Calculating TPM --------"
	#log: "output/counts/tpmcalculator/tpm-calculator.log"
	threads: config["threads"]["TPMCalculator"]
	run:
		# Perform TPM calculation using each sample-level bam file
		shell("TPMCalculator \
			   -k gene_id \
			   -t transcript_id \
			   -o 0 \
			   -g " + config_dict["GTFname"] + " -b {input}")
		# Move output files to proper location
		shell("mv {wildcards.replicate}_genes.* output/counts/tpmcalculator/")

rule mergeTPMCalculator:
	input: expand("output/counts/tpmcalculator/{replicate}_genes.out", replicate=REPLICATE_LIST)
	output: "output/counts/tpmcalculator/tpmcalculator-merged.tsv"
	message: "--------Merging TPMCalculator results--------"
	threads: config["threads"]["mergeTPMCalculator"]
	run:
		shell("./scripts/mergeTPMCalculator.py {input}")

rule featureCounts:
	message: "-----Generating raw counts (featureCounts)-----"
	input: expand("output/{replicate}/bam/{replicate}.bam", replicate=REPLICATE_LIST)
	output: 
		raw = "output/counts/featureCounts/featureCounts.cnt",
		cleaned = "output/counts/featureCounts/featureCount_clean.cnt"
	log: "output/counts/featureCounts/featureCounts.log"
	threads: config["threads"]["featureCount"]
	run:
		# Calculate raw counts from all sample-level bam files
		shell("featureCounts \
				-o {output.raw} \
				-T {threads} \
				-Q 1 \
				-p -M \
				-g gene_id \
				-a " + config_dict["GTFname"] + " {input} \
			2> {log}")
		# Remove the featureCounts header line from output file and re-format the output names
		shell("cat {output.raw} |  egrep -v '#' | \
				sed 's/\Aligned\.sortedByCoord\.out\.bam//g; s/\.bam//g; s/output\/[A-Za-z0-9_-]*\/bam\///g' \
				> {output.cleaned}")

rule HTseq:
	message: "-----Generating raw counts (HTseq)-----"
	input:expand("output/{replicate}/bam/{replicate}.bam", replicate=REPLICATE_LIST)
	output: "output/counts/htseq/htseq-count.tsv"
	log: "output/counts/htseq/HTseq.log"
	threads: config["threads"]["HTseq"]
	run:
		# Compute raw gene-wise counts
		shell("htseq-count \
				 --format bam \
				 --order pos \
				 --mode union \
		 --stranded=no \
				 --type exon \
				 --idattr gene_id \
				 --nprocesses {threads} \
				 --counts_output {output} \
				 {input} " + config_dict["GTFname"] + "  \
				 &> {log}")
		# Add headers to label HTseq tsv output fields
		shell("sed -i '1 i\gene\\t" + "\\t".join(input) + "' {output} &&\
			   sed -i s/'output\/[A-Za-z0-9_-]*\/bam\/'//g {output}")
		# Re-format treatment names in HTseq output
		shell("sed -i 's/\.bam//g' {output}")

rule normalizeFeatureCounts:
	input:
		"output/counts/featureCounts/featureCount_clean.cnt"
	output:
		"output/counts/featureCounts/featureCounts.tpm.tsv",
		"output/counts/featureCounts/featureCounts.fpkm.tsv"
	message: "------- Normalizing featureCounts ---------"
	log: "output/counts/normalize_featureCounts.log"
	threads: config["threads"]["normalizeFeatureCounts"]
	run:
		# Compute TPM and FPKM values featureCounts raw counts
		shell("python ./scripts/normalizeCounts.py featureCounts " + config["GTFname"] + " {input} output/counts/featureCounts/featureCounts " + config["genomeDir"])

rule normalizeHTseq:
	input:
		"output/counts/htseq/htseq-count.tsv"
	output:
		tpm = "output/counts/htseq/htseq-count.tpm.tsv",
		fpkm = "output/counts/htseq/htseq-count.fpkm.tsv"
	message: "------- Normalizing HTseq ---------"
	log: "output/counts/htseq/normalize_HTseq.log"
	threads: config["threads"]["normalizeHTseq"]
	run:
		# Generate exon lengths from gtf file
		shell("./scripts/calc_cdna_len.py " + config["GTFname"] + " gene_id > " + config["genomeDir"] + "cds_length.tsv")
		# Compute TPM and FPKM values HTseq raw counts
		shell("python scripts/normalizeCounts.py HTseq " + config["GTFname"] + " {input} output/counts/htseq/htseq-count " + config["genomeDir"])

rule DEG_featureCounts:
	input: 
		featureCounts = "output/counts/featureCounts/featureCount_clean.cnt"
	output: 
		featureCounts = "output/DEG/featureCount_clean.cnt." + CONTRASTS_FILE.iloc[1][0]  + "_vs_" + CONTRASTS_FILE.iloc[-1][0]  + ".DESeq2.DE_results"
	message: "-------Calculating DEGs {output}---------"
	log: "output/DEG/DEG_featureCounts.log"
	threads: config["threads"]["DEG_featureCounts"]
	run:
		# Compute differentially expressed genes based on contrasts.txt
		shell("run_DE_analysis.pl --matrix {input.featureCounts} --method DESeq2 --samples_file contrasts.txt --output output/DEG")

rule DEG_HTseq:
	input:
		HTseq = "output/counts/htseq/htseq-count.tsv"
	output:
		HTseq = "output/DEG/htseq-count.tsv." + CONTRASTS_FILE.iloc[1][0]  + "_vs_" + CONTRASTS_FILE.iloc[-1][0]  + ".DESeq2.DE_results"
	message: "-------Calculating DEGs {output}---------"
	log: "output/DEG/DEG_HTseq.log"
	threads: config["threads"]["DEG_HTseq"]
	run:
		# Compute differentially expressed genes based on contrasts.txt
		shell("run_DE_analysis.pl --matrix {input.HTseq} --method DESeq2 --samples_file contrasts.txt --output output/DEG")

rule summarize:
	input: 
		tpmcalc = "output/counts/tpmcalculator/tpmcalculator-merged.tsv",
		featureCounts_TPM = "output/counts/featureCounts/featureCounts.tpm.tsv",
		featureCounts_FPKM = "output/counts/featureCounts/featureCounts.fpkm.tsv",
		HTseq_TPM = "output/counts/htseq/htseq-count.tpm.tsv",
		HTseq_FPKM = "output/counts/htseq/htseq-count.fpkm.tsv"
	output: 
		"output/counts/tpmcalculator/tpmcalculator-merged.tsv.average.tsv",
		"output/counts/featureCounts/featureCounts.tpm.tsv.average.tsv",
		"output/counts/htseq/htseq-count.tpm.tsv.average.tsv"
	message: "--------------Summarizing Normalized Counts--------------"
	threads: config["threads"]["summarize"]
	run:
		# Compute summary statistics (average/standard deviation) from TPMcalculator TPM values
		shell("./scripts/summarizeNormalizedCounts.py Gene_Id {input.tpmcalc}")
		# Compute summary statistics (average/standard deviation) from HTseq FPKM/TPM values
		shell("./scripts/summarizeNormalizedCounts.py gene {input.HTseq_TPM}")
		shell("./scripts/summarizeNormalizedCounts.py gene {input.HTseq_FPKM}")
		# Compute summary statistics (average/standard deviation) from featureCounts FPKM/TPM values
		shell("./scripts/summarizeNormalizedCounts.py Geneid {input.featureCounts_TPM}")
		shell("./scripts/summarizeNormalizedCounts.py Geneid {input.featureCounts_FPKM}")
