import sys
import os
import pandas as pd
SAMPLE_FILE = pd.read_table('sraid.txt', sep="\s+", dtype=str).set_index("sample", drop=False)  # enforce str in index
SAMPLE_LIST = SAMPLE_FILE["sample"].values.tolist()

##Create sample output folders
os.makedirs("output/logs/", exist_ok=True)
for path in SAMPLE_LIST:
    os.makedirs("output/" + path + "/stat/", exist_ok=True)
    os.makedirs("output/" + path + "/raw/", exist_ok=True)
    os.makedirs("output/" + path + "/SRA/", exist_ok=True)
	os.makedirs("output/" + path + "/bam/", exist_ok=True)


#REF_FILES = ['TAIR10_chr_all.fas',  'TAIR10_GFF3_genes.gff']
#missing_files = 0
#os.makedirs('ref', exist_ok=True)
#for ref in REF_FILES:
#    if not os.path.exists('ref/' + ref):
#        missing_files = 1
#if missing_files:
#    os.system('curl -L https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff > ref/TAIR10_GFF3_genes.gff; curl -L https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas > ref/TAIR10_chr_all.fas')


configfile: "config.json"

rule all:
    input:
    ### Add step by step output for expand
        expand("output/}sample/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=SAMPLE_LIST)

rule downloadSRA:
    message: "---Downloading SRA files---"
#	input:
		#expand("{reference}", reference=config["reference"])
	threads: config["downloadSRA"]["threads"]
	replicate: config["replicates"]
    output: "output/{sample}/raw/{sample}_{replicate}.fastq.gz"
		#expand("raw_data/{sample}_{replicate}.fastq.gz",
         #   sample=config["samples"],
         #   replicate=config["replicates"])
    log:
        "output/logs/{sample}_fastq-dump.log"
    run:
		shell("prefetch {wildcards.sample} --output-file output/{wildcards.sample}/SRA")
		shell("parallel-fastq-dump --sra-id {wildcards.sample} \
				--threads {threads} --split-3 --gzip \
				--outdir output/{wildcards.sample}/raw \
				2> {log}"
			  
		
		
        #for sample in config["samples"]:
        #    shell("prefetch " + sample)
        #    shell("parallel-fastq-dump --sra-id " + sample +
        #          " --threads " + str(config["maxThreads"]) +
        #          " --split-3 --outdir ./raw_data  --gzip 2> {log}")

rule trim:
    input:
        expand("raw_data/{sample}_{replicate}.fastq.gz",
            sample=config["samples"],
            replicate=config["replicates"])
    output:
        expand("trim/{sample}_{replicate}.fq.gz",
        sample=config["samples"],
        replicate=config["replicates"])
    log:
        "logs/trim.log"
    threads:
        config["maxThreads"]
    run:
        for sample in config["samples"]:
                shell("trim_galore --paired --three_prime_clip_R1 5 --three_prime_clip_R2" +
                      " 5 --cores 2 --max_n 40 --gzip -o trim raw_data/" + sample +
                      "_1.fastq.gz raw_data/" + sample + "_2.fastq.gz 2> {log}")
        shell("for file in trim/*_val_*; do mv $file $(echo $file | sed s/_val_[0-9]//); done")
        shell("fastqc " + " ".join(expand("trim/{sample}_{replicate}.fq.gz",sample=config["samples"],replicate=config["replicates"])))

rule sample:
    input:
        expand("trim/{sample}_{replicate}.fq.gz",
                sample=config["samples"],
                replicate=config["replicates"])
    output:
        expand("trim/sample/{sample}_{replicate}_1k.fq.gz",
            sample=config["samples"],
            replicate=config["replicates"])
    run:
        for sample in config["samples"]:
            shell("zcat < trim/" + sample + "_1.fq.gz | seqkit sample -n 1000 -o trim/sample/" + sample + "_1_1k.fq.gz")
            shell("zcat < trim/" + sample + "_2.fq.gz | seqkit sample -n 1000 -o trim/sample/" + sample + "_2_1k.fq.gz")

rule align:
    input:
        expand("trim/sample/{sample}_{replicate}_1k.fq.gz",
            sample=config["samples"],
            replicate=config["replicates"])
    output:
        expand("bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=config["samples"])
    log:
        "logs/align.log"
    threads:
        config["maxThreads"]
    run:
        for sample in config["samples"]:
            shell("STAR --runMode alignReads --runThreadN 1" +
                  " --readFilesCommand gunzip -c --outFilterMultimapNmax 10 --alignIntronMin 25" +
                  " --alignIntronMax 10000 --genomeDir " + config["genomeDir"] +
                  " --readFilesIn trim/sample/" + sample + "_1_1k.fq.gz trim/sample/" + sample + "_2_1k.fq.gz" +
                  " --outSAMtype BAM SortedByCoordinate --outFileNamePrefix " + config["bamPath"] + sample + ".bam")
