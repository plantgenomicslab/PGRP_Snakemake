import os
import sys

import pandas as pd

# Load run and sample information from the sraRunsByExperiment.tsv input file (user must provide)
SAMPLES_FILE = pd.read_csv("RunsByExperiment.tsv", sep="\t")
SAMPLE_LIST = list(set(SAMPLES_FILE["Run"].values.tolist()))
REPLICATE_LIST = list(set(SAMPLES_FILE["Replicate"].tolist()))
REPLICATE_LOOKUP = SAMPLES_FILE.groupby("Replicate")['Run'].apply(list).to_dict()

# Load configs (Snakemake's configfile: directive populates the `config` dict;
# no need to re-load yaml manually).
configfile: "config.yml"

SOURCE = config["source"]
LAYOUT = config["libraryType"] #TODO enable mixed library types
if LAYOUT == "PAIRED":
	PAIR_LIST = config["PAIR_LIST"]
else:
	PAIR_LIST = ["_1", "_2"]

# BBDuk pre-alignment contaminant filter (optional). Enabled when config["bbduk_enable"] = true.
# References must be built once via:  bash scripts/build_bbduk_refs.sh
BBDUK_ENABLE = config.get("bbduk_enable", False)
if BBDUK_ENABLE:
	include: "rules/bbduk_filter.smk"

def _trim_inputs(wc):
	if LAYOUT == "PAIRED":
		return {
			"fwd_fastq": f"output/{wc.replicate}/{wc.sample}/trim/{wc.sample}{PAIR_LIST[0]}_trimmed.fq.gz",
			"rev_fastq": f"output/{wc.replicate}/{wc.sample}/trim/{wc.sample}{PAIR_LIST[1]}_trimmed.fq.gz",
		}
	return {"fwd_fastq": f"output/{wc.replicate}/{wc.sample}/trim/{wc.sample}_trimmed.fq.gz"}

def _bbduk_inputs(wc):
	if LAYOUT == "PAIRED":
		return {
			"fwd_fastq": f"output/{wc.replicate}/{wc.sample}/bbduk/{wc.sample}{PAIR_LIST[0]}_clean.fq.gz",
			"rev_fastq": f"output/{wc.replicate}/{wc.sample}/bbduk/{wc.sample}{PAIR_LIST[1]}_clean.fq.gz",
		}
	return {"fwd_fastq": f"output/{wc.replicate}/{wc.sample}/bbduk/{wc.sample}_clean.fq.gz"}

def align_inputs(wc):
	return _bbduk_inputs(wc) if BBDUK_ENABLE else _trim_inputs(wc)

# Check for an indexed reference genome or prepared genome for RSEM
if "RSEM" in config["readCounting"]:
	if not os.path.exists(f"{config['RSEM_prepared_genome']}.seq"):
		sys.exit(f"Cannot locate '{config['RSEM_prepared_genome']}'.\nProvide a reference genome prepared with 'rsem-prepare-reference' to use RSEM.\nExiting...")

for ref in config["ref"]:
	if not os.path.exists(f"{config['genomeDir']}/{ref}"):
		sys.exit("The " + ref + " file is missing.\nHave you provided the correct paths to a reference genome indexed by STAR?\nExiting...")

# Configure DEG calculation (user must provide)
#TODO: write better method for identifying DEG contrasts output files. Maybe 'checkpoint'?
contrasts = os.path.exists(config["sample_contrast"])
if config["runDEG"]:
	try:
		CONTRASTS_FILE = pd.read_csv(config["sample_contrast"], sep="\t", header=None)
		cTop = CONTRASTS_FILE.iloc[-1][0]
		cBottom = CONTRASTS_FILE.iloc[-1][1]
	except (FileNotFoundError, pd.errors.EmptyDataError):
		cTop = ""
		cBottom = ""
else:
	CONTRASTS_FILE = pd.DataFrame([["",""],["",""]])
	cTop = ""
	cBottom = ""

# create the Sample tab-delimited text file indicating biological replicate relationships
replicate_relationship= ""
for index, row in SAMPLES_FILE.iterrows():
	replicate_relationship += f"{row['Treatment']}\t{row['Replicate']}\n"
	with open('replication_relationship.txt', 'w') as f:
		f.write(replicate_relationship)

# Create sample output folders
os.makedirs("output/logs/", exist_ok=True)
os.makedirs("output/DEG/", exist_ok=True)

counting_options  = ("RSEM","featureCounts","HTseq","TPMcalculator")
for opt in counting_options:
	if opt in config["readCounting"]:
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
	rsem = False
	star = False
	
	# Count and DEG files
	if "RSEM" in config["readCounting"]:
		inputs += ["output/counts/RSEM/RSEM_TPM.tsv.average.tsv",
			   "output/counts/RSEM/RSEM_TPM.tsv"]
		if config["runDEG"]:
			inputs += ["output/DEG/RSEM_expected_count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results",
				   "output/DEG/RSEM_expected_count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset"]
	if "featureCounts" in config["readCounting"]:
		inputs += ["output/counts/featureCounts/featureCounts.cnt",
			   "output/counts/featureCounts/featureCounts.tpm.tsv",
			   "output/counts/featureCounts/featureCounts.fpkm.tsv",
			   "output/counts/featureCounts/featureCounts.tpm.tsv.average.tsv",
			   "output/counts/featureCounts/featureCount_clean.cnt"]
		if config["runDEG"]:
			inputs += ["output/DEG/featureCount_clean.cnt." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset",
				   "output/DEG/featureCount_clean.cnt." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results"]
	if "HTseq" in config["readCounting"]:
		inputs += ["output/counts/htseq/htseq-count.tsv",
			   "output/counts/htseq/htseq-count.tpm.tsv",
			   "output/counts/htseq/htseq-count.fpkm.tsv",
			   "output/counts/htseq/htseq-count.tpm.tsv.average.tsv"]
		if config["runDEG"]:
			inputs += ["output/DEG/htseq-count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results",
				   "output/DEG/htseq-count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset"]
	if "TPMcalculator" in config["readCounting"]:
		inputs += ["output/counts/tpmcalculator/tpmcalculator-merged.tsv",
			   "output/counts/tpmcalculator/tpmcalculator-merged.tsv.average.tsv"]
		
	# Per-replicate data and processing files
	for replicate in REPLICATE_LIST:	
		if "TPMcalculator" in config["readCounting"]:
			inputs.append("output/counts/tpmcalculator/" + replicate + "_genes.out")
			star = True

		if "RSEM" in config["readCounting"]:
			inputs.append("output/counts/RSEM/" + replicate + ".genes.results")
			inputs.append("output/" + replicate + "/bam/" + replicate + ".RSEM.bam")
			rsem = True
		
		if ("featureCounts" in config["readCounting"]) | ("HTseq" in config["readCounting"]) | ("TPMcalculator" in config["readCounting"]):
				inputs.append("output/" + replicate + "/bam/" + replicate + ".STAR.bam")
				star = True

		for sample in REPLICATE_LOOKUP[replicate]:
			if star:
				inputs.append("output/" + replicate + "/" + sample + "/bam/" + sample + ".bamAligned.sortedByCoord.out.bam")
			if rsem:
				inputs.append("output/" + replicate + "/" + sample + "/bam/" + sample + ".xs.bamAligned.toTranscriptome.out.bam")

			if LAYOUT == "PAIRED":
				for pair in PAIR_LIST:
					inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample  + pair + ".fastq.gz")
					inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample  + pair + "_fastqc.zip")
					inputs.append("output/" + replicate + "/" + sample + "/trim/" + sample + pair + "_trimmed.fq.gz")
			elif LAYOUT == "SINGLE":
				inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample + ".fastq.gz")
				inputs.append("output/" + replicate + "/" + sample + "/raw/" + sample + "_fastqc.zip")
				inputs.append("output/" + replicate + "/" + sample + "/trim/" + sample + "_trimmed.fq.gz")
			else:
				print("ERROR: libraryType must be 'SINGLE' or 'PAIRED' ... see config/config.json to set libraryType")
				quit()
	return(inputs)

if LAYOUT == "PAIRED":
	ruleorder: trim_PAIRED > trim_SINGLE
	ruleorder: convertSRAtoFastq_PAIRED > convertSRAtoFastq_SINGLE
	ruleorder: importRaw_PAIRED > importRaw_SINGLE
	ruleorder: fastqc_raw_PAIRED > fastqc_raw_SINGLE
	if SOURCE == "local":
		ruleorder: importRaw_PAIRED > convertSRAtoFastq_PAIRED
	else:
		ruleorder: convertSRAtoFastq_PAIRED > importRaw_PAIRED
elif LAYOUT == "SINGLE":
	ruleorder: trim_SINGLE > trim_PAIRED
	ruleorder: convertSRAtoFastq_SINGLE > convertSRAtoFastq_PAIRED
	ruleorder: importRaw_SINGLE > importRaw_PAIRED
	ruleorder: fastqc_raw_SINGLE > fastqc_raw_PAIRED
	if SOURCE == "local":
		ruleorder: importRaw_SINGLE > convertSRAtoFastq_SINGLE
	else:
		ruleorder: convertSRAtoFastq_SINGLE > importRaw_SINGLE

localrules: all, importRaw_PAIRED, importRaw_SINGLE


def flagstat_check(output, log):
	"""Run `samtools flagstat` against a BAM and abort the rule on failure."""
	shell(
		"echo '--------Checking " + str(output) + "----------' && "
		"samtools flagstat " + str(output) + " "
		"|| { echo 'samtools flagstat found errors in " + str(output) + ". "
		"Check log here: " + str(log) + ". Exiting......' >&2 ; exit 1; }"
	)

rule all:
	input:
		allInput()


# -------------------------------------------------------------------------
# Rule modules (split out for readability; each file holds a thematic group
# of rules. Helpers used across modules — flagstat_check, allInput,
# align_inputs — live above this point in this file).
# -------------------------------------------------------------------------
include: "rules/import.smk"
include: "rules/qc.smk"
include: "rules/align.smk"
include: "rules/count.smk"
include: "rules/deg.smk"
include: "rules/summarize.smk"
