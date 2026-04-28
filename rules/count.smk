rule TPMCalculator:
	input: "output/{replicate}/bam/{replicate}.STAR.bam"
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
			   -g " + config["GTFname"] + " -b {input}")
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
	input: expand("output/{replicate}/bam/{replicate}.STAR.bam", replicate=REPLICATE_LIST)
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
				-a " + config["GTFname"] + " {input} \
			2> {log}")
		# Remove the featureCounts header line from output file and re-format the output names
		shell("cat {output.raw} |  egrep -v '#' | \
				sed 's/\Aligned\.sortedByCoord\.out\.bam//g; s/\.bam//g; s/output\/[A-Za-z0-9_-]*\/bam\///g' \
				> {output.cleaned}")

rule HTseq:
	message: "-----Generating raw counts (HTseq)-----"
	input:expand("output/{replicate}/bam/{replicate}.STAR.bam", replicate=REPLICATE_LIST)
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
				 {input} " + config["GTFname"] + "  \
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

