rule importRaw_PAIRED:
	output: 
		one = "output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[0] + ".fastq.gz",
		two = "output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[1] + ".fastq.gz"
	message: "Importing raw data: {wildcards.sample}"
	run:
		shell("ln -s " + config["rawInputDir"] + "/{wildcards.sample}" + PAIR_LIST[0] + ".fastq.gz {output.one}")
		shell("ln -s " + config["rawInputDir"] + "/{wildcards.sample}" + PAIR_LIST[1] + ".fastq.gz {output.two}")

rule importRaw_SINGLE:
	output: "output/{replicate}/{sample}/raw/{sample}.fastq.gz"
	message: "Importing raw data: {wildcards.sample}"
	run:
		shell("ln -s " + config["rawInputDir"] + "/{wildcards.sample}.fastq.gz {output}")

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
		"output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[0] + ".fastq.gz",
		"output/{replicate}/{sample}/raw/{sample}" + PAIR_LIST[1] + ".fastq.gz"
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

