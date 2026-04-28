rule summarizeTPMCalc:
	input: tpmcalc = "output/counts/tpmcalculator/tpmcalculator-merged.tsv"
	output: "output/counts/tpmcalculator/tpmcalculator-merged.tsv.average.tsv"
	message: "--------------Summarizing Normalized Counts--------------"
	run:
		# Compute summary statistics (average/standard deviation) from TPMcalculator TPM values
		shell("./scripts/summarizeNormalizedCounts.py Gene_Id {input.tpmcalc}")

rule summarizeHTseq:
	input:
		HTseq_TPM = "output/counts/htseq/htseq-count.tpm.tsv",
		HTseq_FPKM = "output/counts/htseq/htseq-count.fpkm.tsv"
	output: "output/counts/htseq/htseq-count.tpm.tsv.average.tsv"
	message: "--------------Summarizing HTseq--------------"
	run:
		# Compute summary statistics (average/standard deviation) from HTseq FPKM/TPM values
		shell("./scripts/summarizeNormalizedCounts.py gene {input.HTseq_TPM}")
		shell("./scripts/summarizeNormalizedCounts.py gene {input.HTseq_FPKM}")

rule summarizeFeatureCounts:
	input:
		featureCounts_TPM = "output/counts/featureCounts/featureCounts.tpm.tsv",
		featureCounts_FPKM = "output/counts/featureCounts/featureCounts.fpkm.tsv"
	output: "output/counts/featureCounts/featureCounts.tpm.tsv.average.tsv"
	message: "--------------Summarizing featureCounts--------------"
	run:
		# Compute summary statistics (average/standard deviation) from featureCounts FPKM/TPM values
		shell("./scripts/summarizeNormalizedCounts.py Geneid {input.featureCounts_TPM}")
		shell("./scripts/summarizeNormalizedCounts.py Geneid {input.featureCounts_FPKM}")

rule summarizeRSEM:
	input: input_DEG_RSEM
	output: 
		"output/counts/RSEM/RSEM_TPM.tsv.average.tsv",
		"output/counts/RSEM/RSEM_TPM.tsv"
	message: "------ Summarizing RSEM ------"
	run:
		# Make RSEM tpm/fpkm matrices
		shell("python ./scripts/makeRSEMMatrix.py RunsByExperiment.tsv output/counts/RSEM TPM")
		shell("python ./scripts/makeRSEMMatrix.py RunsByExperiment.tsv output/counts/RSEM FPKM")
		# Summarize RSEM data
		shell("./scripts/summarizeNormalizedCounts.py gene_id output/counts/RSEM/RSEM_TPM.tsv")
		shell("./scripts/summarizeNormalizedCounts.py gene_id output/counts/RSEM/RSEM_FPKM.tsv")
