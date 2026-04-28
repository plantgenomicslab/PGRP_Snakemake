rule align:
	input:
		unpack(align_inputs)
	output: "output/{replicate}/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam"
	message: "-----Aligning {wildcards.sample}-----"
	log: "output/{replicate}/{sample}/logs/{sample}_align.log"
	threads: config["threads"]["align"]
	run:
		# {input} expands to either fwd_fastq or fwd_fastq+rev_fastq depending on LAYOUT (see align_inputs).
		shell("STAR --runMode alignReads --runThreadN {threads} \
				--outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 \
				--quantMode TranscriptomeSAM GeneCounts \
				--outBAMsortingBinsN 200 \
				--genomeDir " + config["genomeDir"]  + " \
				--readFilesCommand gunzip -c --readFilesIn {input} \
				--outSAMtype BAM SortedByCoordinate --outFileNamePrefix output/{wildcards.replicate}/{wildcards.sample}/bam/{wildcards.sample}.bam  \
				2> {log}")
		flagstat_check(output, log)

rule alignRSEM:
	input:
		unpack(align_inputs)
	output: "output/{replicate}/{sample}/bam/{sample}.xs.bamAligned.toTranscriptome.out.bam"
	message: "-----Aligning for RSEM: {wildcards.sample}-----"
	log: "output/{replicate}/{sample}/logs/{sample}_alignRSEM.log"
	threads: config["threads"]["align"]
	run:
		shell("STAR --runMode alignReads --runThreadN {threads} \
			--outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 25000 \
			--limitBAMsortRAM 20000000000 \
			--outBAMsortingBinsN 200 \
			--genomeDir " + config["genomeDir"]  + " \
			--readFilesCommand gunzip -c --readFilesIn {input} \
			--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
			--quantTranscriptomeBan IndelSoftclipSingleend  \
			--alignEndsType EndToEnd  \
			--outFileNamePrefix output/{wildcards.replicate}/{wildcards.sample}/bam/{wildcards.sample}.xs.bam")

def mergeInputRSEM(wildcards):
	inputs = expand("output/" + wildcards.replicate + "/{sample}/bam/{sample}.xs.bamAligned.toTranscriptome.out.bam", sample=REPLICATE_LOOKUP[wildcards.replicate])
	return(inputs)

rule mergeRSEM:
	input: mergeInputRSEM
	output: "output/{replicate}/bam/{replicate}.RSEM.bam"
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
			shell("ln -s --relative {input} {output}")
			print("1 run per sample...skipping BAM merge")
		# Perform check on output bam file to ensure it is not corrupted
		flagstat_check(output, log)

#TODO: Need to merge bam files from each run to calculate replicate-level expression
rule calculateRSEMExpression:
	input: "output/{replicate}/bam/{replicate}.RSEM.bam"
	output: "output/counts/RSEM/{replicate}.genes.results"
	message: "-----  Calculating RSEM expression: {wildcards.replicate}------"
	threads: config["threads"]["calculateRSEMExpression"]
	params:
		paired = "--paired-end" if LAYOUT == "PAIRED" else "",
	run:
		shell("rsem-calculate-expression --bam --no-bam-output -p {threads} \
			{params.paired} {input} " + config["RSEM_prepared_genome"] + " output/counts/RSEM/{wildcards.replicate}")

def mergeInputSTAR(wildcards):
	inputs = expand("output/" + wildcards.replicate + "/{sample}/bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=REPLICATE_LOOKUP[wildcards.replicate])
	return(inputs)

rule mergeSTAR:
	input: mergeInputSTAR
	output: "output/{replicate}/bam/{replicate}.STAR.bam"
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
			shell("ln -s {input} {output}")
			print("1 run per sample...skipping BAM merge")
		# Perform check on output bam file to ensure it is not corrupted
		flagstat_check(output, log)

