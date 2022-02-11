import sys
import os
import pandas as pd

SAMPLE_FILE = pd.read_table("sraid.txt", sep="\s+", dtype=str).set_index("sample", drop=False)  # enforce str in index
SAMPLE_LIST = SAMPLE_FILE["sample"].values.tolist()

REPLICATE_LIST=["1", "2"]

TREATMENT_FILE = pd.read_csv("sraRunsbyExperiment.csv")
TREATMENT_LIST = list(set(TREATMENT_FILE["Treatment"].values.tolist()))
TREATMENT_LOOKUP = TREATMENT_FILE.groupby("Treatment")['Run'].apply(list).to_dict()

# Create sample output folders
os.makedirs("output/counts/", exist_ok=True)
os.makedirs("output/sra/", exist_ok=True)
os.makedirs("output/sra/logs/", exist_ok=True)
os.makedirs("output/logs/", exist_ok=True)

for path in SAMPLE_LIST:
    os.makedirs("output/" + path + "/raw/", exist_ok=True)
    os.makedirs("output/" + path + "/bam/", exist_ok=True)
    os.makedirs("output/" + path + "/logs/", exist_ok=True)
    os.makedirs("output/" + path + "/trim/", exist_ok=True)

# Load the config file as a dictionary
cf = open("config.json")
config_dict = json.load(cf)
cf.close()

# Check for an indexed reference genome based on required files in config.json
for ref in config_dict["ref"]:
    if not os.path.exists(ref):
        print("The " + ref + " file is missing. Have you downloaded and indexed a refernce genome? Terminating...")
        exit()

configfile: "config.json"

def allInput():
    inputs = ["output/counts/featureCounts.cnt",
              "output/counts/htseq-count.tsv"]
    for treatment in TREATMENT_LIST:
        inputs.append("output/" + treatment + "/bam/" + treatment + ".bam")
        for sample in TREATMENT_LOOKUP[treatment]:
            inputs.append("output/" + treatment + "/" + sample + "/bam/" + sample + ".bamAligned.sortedByCoord.out.bam")
            for replicate in REPLICATE_LIST:
                inputs.append("output/" + treatment + "/" + sample + "/raw/" + sample + "_" + replicate + ".fastq.gz")
                inputs.append("output/" + treatment + "/" + sample + "/raw/" + sample + "_" + replicate + "_fastqc.zip")
                inputs.append("output/" + treatment + "/" + sample + "/trim/" + sample + "_" + replicate + ".fq.gz")
    return(inputs)

rule all:
    input:
        allInput()

rule fetchSRA:
    output: "output/sra/{sample}.sra"
    message: "-----Fetching {wildcards.sample} SRA files-----"
    threads: config["threads"]["fetchSRA"]
    log: "output/sra/logs/{sample}_downloadSRA.log"
    run:
        shell("prefetch {wildcards.sample} --output-file {output} 2> {log}")
        # Check sra files to make sure they are valid
        shell("echo '--------Validating {wildcards.sample}.sra--------'")
        shell("set +e")
        shell("if ! vdb-validate {output}; then echo 'vdb-validate found errors in {output}. Check log here:{log}. Exiting......' && exit 1; fi")

rule convertSRAtoFastq:
    input: "output/sra/{sample}.sra"
    output:
        "output/{treatment}/{sample}/raw/{sample}_1.fastq.gz",
        "output/{treatment}/{sample}/raw/{sample}_2.fastq.gz"
    message: "-----Converting {wildcards.sample} SRA to Fastq files-----"
    threads: config["threads"]["convertSRAtoFastq"]
    log: "output/{treatment}/{sample}/logs/{sample}_fastqdump.log"
    run:
        shell("parallel-fastq-dump --sra-id {wildcards.sample} \
				--threads {threads} --split-e --gzip \
				--outdir output/{wildcards.treatment}/{wildcards.sample}/raw \
				2> {log}")

rule fastqc_raw:
    input:
        "output/{treatment}/{sample}/raw/{sample}_1.fastq.gz",
        "output/{treatment}/{sample}/raw/{sample}_2.fastq.gz"
    output:
        "output/{treatment}/{sample}/raw/{sample}_1_fastqc.zip",
        "output/{treatment}/{sample}/raw/{sample}_2_fastqc.zip"
    message: "-----Running Fastqc_raw {wildcards.sample}-----"
    threads: config["threads"]["fastqc_raw"]
    log: "output/{treatment}/{sample}/logs/{sample}_raw_fastqc.log"
    run:
        shell("fastqc --threads {threads} {input} 2> {log}")

rule trim:
    input:
        fwd_fastq = "output/{treatment}/{sample}/raw/{sample}_1.fastq.gz",
        rev_fastq = "output/{treatment}/{sample}/raw/{sample}_2.fastq.gz"
    output:
        "output/{treatment}/{sample}/trim/{sample}_1.fq.gz",
        "output/{treatment}/{sample}/trim/{sample}_2.fq.gz"
    message: "-----Trimming {wildcards.sample}-----"
    log: "output/{treatment}/{sample}/logs/{sample}_trim.log"
    threads: config["threads"]["trim"]
    run:
        shell("trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
             --cores {threads} --max_n 40 --gzip -o output/{wildcards.treatment}/{wildcards.sample}/trim {input.fwd_fastq} {input.rev_fastq} \
             2> {log}")
        shell("for file in output/{wildcards.treatment}/{wildcards.sample}/trim/*_val_*; do mv $file $(echo $file | sed s/_val_[0-9]//); done")
        shell("fastqc --threads {threads} {output}")

rule align:
    input:
        fwd_fastq = "output/{treatment}/{sample}/trim/{sample}_1.fq.gz",
        rev_fastq = "output/{treatment}/{sample}/trim/{sample}_2.fq.gz"
    output:"output/{treatment}/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam"
    message: "-----Aligning {wildcards.sample}-----"
    log: "output/{treatment}/{sample}/logs/{sample}_align.log"
    threads: config["threads"]["align"]
    run:
        shell("STAR --runMode alignReads --runThreadN {threads} \
                --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 10000 \
		--genomeDir " + config["genomeDir"]  + " \
	        --readFilesCommand gunzip -c --readFilesIn {input.fwd_fastq} {input.rev_fastq} \
                --outSAMtype BAM SortedByCoordinate --outFileNamePrefix output/{treatment}/{wildcards.sample}/bam/{wildcards.sample}.bam  \
		2> {log}")
        # Perform quick check on output bam file to ensure it is not corrupted
        shell("echo '--------Checking {output}----------'")
        shell("set +e")
        shell("if ! samtools quickcheck {output}; then echo 'samtools quickcheck found errors in {output}. Check log here: {log}. Exiting......' && exit 1; fi")

def mergeInput(wildcards):
    inputs = expand("output/" + wildcards.treatment + "/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=TREATMENT_LOOKUP[wildcards.treatment])
    return(inputs)

rule merge:
    input: mergeInput
    output: "output/{treatment}/bam/{treatment}.bam"
    message: "-----Merging bam files by experiment-----"
    log: "output/{treatment}/bam/{treatment}_merge.log"
    threads: config["threads"]["merge"]
    run:
        shell("sambamba merge --nthreads {threads} {output} {input}")
        # Create bam index folder for HTseq
        shell("samtools index {output}")

rule featureCounts:
    message: "-----Generating raw counts (featureCounts)-----"
    input: expand("output/{treatment}/bam/{treatment}.bam", treatment=TREATMENT_LIST)
    output: "output/counts/featureCounts.cnt"
    log: "output/counts/featureCounts.log"
    threads: config["threads"]["featureCount"]
    run:
        shell("featureCounts -o output/counts/featureCounts.cnt  -T {threads} -p -a " + config_dict["GTFname"] + " {input} \
    		2> {log}")

rule HTseq:
    message: "-----Generating raw counts (HTseq)-----"
    input:expand("output/{treatment}/bam/{treatment}.bam", treatment=TREATMENT_LIST)
    output: "output/counts/htseq-count.tsv"
    log: "output/counts/HTseq.log"
    threads: config["threads"]["HTseq"]
    run:
        # Compute raw gene-wise counts
        shell("htseq-count \
                 --format bam \
                 --order pos \
                 --mode union \
                 --nprocesses {threads} \
                 --counts_output {output} \
                 {input} " + config_dict["GTFname"] + "  \
                 2> {log}")
        # Add headers to label HTseq tsv output fields
        shell("sed -i '1 i\gene\\t" + "\\t".join(input) + "' output/counts/htseq-count.tsv &&\
               sed -i s/'output\/[A-Za-z0-9_]*\/bam\/'//g output/counts/htseq-count.tsv")

rule normalizeFeatureCounts:
    input:
        "output/counts/featureCounts.cnt"
    output:
        "output/counts/featureCounts.tpm.tsv",
        "output/counts/featureCounts.fpkm.tsv"
    message: "------- Normalizing featureCounts ---------"
    log: "output/counts/normalize_featureCounts.log"
    threads: config["threads"]["normalizeFeatureCounts"]
    run:
        shell("cat {input} |  egrep -v '#' | \
               sed 's/\Aligned\.sortedByCoord\.out\.bam//g; s/\.bam//g' \
               > output/counts/featureCount_clean.cnt")
        shell("./normalizeCounts.py featureCounts " + config["GTFname"] + " output/counts/featureCount_clean.cnt output/counts/featureCounts")

rule normalizeHTseq:
    input:
        "output/counts/htseq-count.tsv"
    output:
        "output/counts/htseq-count.tpm.tsv",
        "output/counts/htseq-count.fpkm.tsv"
    message: "------- Normalizing HTseq ---------"
    log: "output/counts/normalize_HTseq.log"
    threads: config["threads"]["normalizeHTseq"]
    run:
        shell("sed {input} 's/\.bam//g' > output/counts/htseq-count_clean.cnt")
        shell("./normalizeCounts.py HTseq " + config["GTFname"] + " output/counts/htseq-count_clean.cnt output/counts/htseq-count")
