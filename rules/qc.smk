rule fastqc_raw_PAIRED:
	input:
		"output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[0] + ".fastq.gz",
		"output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[1] + ".fastq.gz"
	output:
		"output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[0] + "_fastqc.zip",
		"output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[1] + "_fastqc.zip"
	message: "-----Running Fastqc_raw {wildcards.sample}-----"
	threads: config["threads"]["fastqc_raw"]
	log: "output/{replicate}/{sample}/logs/{sample}_raw_fastqc.log"
	run:
		# Run quality control on raw reads
		shell("fastqc --threads {threads} {input} --outdir output/{wildcards.replicate}/{wildcards.sample}/raw/ 2> {log}")

rule fastqc_raw_SINGLE:
	input: "output/{replicate}/{sample}/raw/{sample}.fastq.gz"
	output: "output/{replicate}/{sample}/raw/{sample}_fastqc.zip"
	message: "-----Running Fastqc_raw {wildcards.sample}-----"
	threads: config["threads"]["fastqc_raw"]
	log: "output/{replicate}/{sample}/logs/{sample}_raw_fastqc.log"
	run:
		# Run quality control on raw reads
		shell("fastqc --threads {threads} {input} --outdir output/{wildcards.replicate}/{wildcards.sample}/raw/ 2> {log}")

rule trim_PAIRED:
	input:
		fwd_fastq = "output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[0] + ".fastq.gz",
		rev_fastq = "output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[1] + ".fastq.gz"
	output:
		fwd_fastq = "output/{replicate}/{sample}/trim/{sample}" + PAIR_LIST[0] + "_trimmed.fq.gz",
		rev_fastq = "output/{replicate}/{sample}/trim/{sample}" + PAIR_LIST[1] + "_trimmed.fq.gz"
	message: "-----Trimming {wildcards.sample}-----"
	log: "output/{replicate}/{sample}/logs/{sample}_trim.log"
	threads: config["threads"]["trim"]
	run:
		# Trim reads using fastp
		shell("fastp --detect_adapter_for_pe \
		--overrepresentation_analysis \
		--cut_right \
		--thread {threads} \
		--html output/{wildcards.replicate}/{wildcards.sample}/trim/{wildcards.sample}.fastp.html \
		--json output/{wildcards.replicate}/{wildcards.sample}/trim/{wildcards.sample}.fastp.json \
		-i {input.fwd_fastq} -I  {input.rev_fastq} \
		-o {output.fwd_fastq} -O {output.rev_fastq} 2> {log}")
		# Run quality control on trimmed reads
		shell("fastqc --threads {threads} {output} --outdir output/{wildcards.replicate}/{wildcards.sample}/trim 2> {log}")

rule trim_SINGLE:
	input: "output/{replicate}/{sample}/raw/{sample}.fastq.gz"
	output:"output/{replicate}/{sample}/trim/{sample}_trimmed.fq.gz"
	message: "-----Trimming {wildcards.sample}-----"
	log: "output/{replicate}/{sample}/logs/{sample}_trim.log"
	threads: config["threads"]["trim"]
	run:
		# Trim reads using fastp
		shell("fastp  \
		--overrepresentation_analysis \
		--cut_right \
		--thread {threads} \
		--html output/{wildcards.replicate}/{wildcards.sample}/trim/{wildcards.sample}.fastp.html \
		--json output/{wildcards.replicate}/{wildcards.sample}/trim/{wildcards.sample}.fastp.json \
		-i {input} \
		-o {output} 2> {log}")
		# Run quality control on trimmed reads
		shell("fastqc --threads {threads} {output} --outdir output/{wildcards.replicate}/{wildcards.sample}/trim 2> {log}")

