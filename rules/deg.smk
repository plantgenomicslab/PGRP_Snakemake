rule DEG_featureCounts:
	input: 
		featureCounts = "output/counts/featureCounts/featureCount_clean.cnt"
	output: 
		featureCounts = "output/DEG/featureCount_clean.cnt." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results",
		featureCounts_subset = "output/DEG/featureCount_clean.cnt." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset",
		PtR = "output/counts/featureCounts/featureCount_clean.cnt.minRow10.CPM.log2.centered.prcomp.principal_components.pdf"
	message: "-------Calculating DEGs {output}---------"
	log: "output/DEG/DEG_featureCounts.log"
	threads: config["threads"]["DEG_featureCounts"]
	params:
		replication = os.path.join(os.getcwd(), config["rep_relations"]),
		matrix =  os.path.join(os.getcwd(), "output/counts/featureCounts/featureCount_clean.cnt")
	run:
		# Compute differentially expressed genes based on deg_samples.txt
		shell("run_DE_analysis.pl --matrix {input.featureCounts} --method DESeq2 --samples_file " + config["rep_relations"] + " --contrasts " + config["sample_contrast"] + " --output output/DEG")
		shell("cd output/DEG && analyze_diff_expr.pl --samples {params.replication} --matrix {params.matrix} -P 0.001 -C 2")
		shell("cd output/DEG && analyze_diff_expr.pl --samples {params.replication} --matrix {params.matrix} -P 0.01 -C 1")
		shell("cd output/counts/featureCounts && PtR --matrix featureCount_clean.cnt --min_rowSums 10 -s {params.replication}  --log2 --CPM --sample_cor_matrix --CPM --center_rows --prin_comp 3")


rule DEG_HTseq:
	input:
		HTseq = "output/counts/htseq/htseq-count.tsv"
	output:
		HTseq = "output/DEG/htseq-count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results",
		HTSeq_subset =  "output/DEG/htseq-count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset",
		PtR = "output/counts/htseq/htseq-count.tsv.minRow10.CPM.log2.centered.prcomp.principal_components.pdf"
	message: "-------Calculating DEGs {output}---------"
	log: "output/DEG/DEG_HTseq.log"
	params:
		 replication = os.path.join(os.getcwd(), config["rep_relations"]),
		 matrix =  os.path.join(os.getcwd(), "output/counts/htseq/htseq-count.tsv")
	threads: config["threads"]["DEG_HTseq"]
	run:
		# Compute differentially expressed genes based on deg_samples.txt
		shell("run_DE_analysis.pl --matrix {input.HTseq} --method DESeq2 --samples_file" + config["rep_relations"] + "--contrasts " + config["sample_contrast"] + " --output output/DEG")
		shell("cd output/DEG && analyze_diff_expr.pl --samples {params.replication} --matrix {params.matrix} -P 0.001 -C 2")
		shell("cd output/DEG && analyze_diff_expr.pl --samples {params.replication} --matrix {params.matrix} -P 0.01 -C 1")
		shell("cd output/counts/htseq && PtR --matrix htseq-count.tsv --min_rowSums 10 -s {params.replication}  --log2 --CPM --sample_cor_matrix --CPM --center_rows --prin_comp 3")

def input_DEG_RSEM(wildcards):
	inputs = expand("output/counts/RSEM/{replicate}.genes.results", replicate=REPLICATE_LIST)
	return(inputs)

rule DEG_RSEM:
	input: input_DEG_RSEM
	output: "output/DEG/RSEM_expected_count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results",
		"output/DEG/RSEM_expected_count.tsv." + cTop  + "_vs_" + cBottom  + ".DESeq2.DE_results.P0.01_C1.DE.subset",
		"output/counts/RSEM/RSEM_expected_count.tsv.minRow10.CPM.log2.centered.prcomp.principal_components.pdf"
	message: "-------Calculating RSEM DEGs {output}---------"
	log: "output/DEG/DEG_RSEM.log"
	threads: config["threads"]["DEG_RSEM"]
	params:
		rep_relations = os.path.join(os.getcwd(), config["rep_relations"]),
		matrix =  os.path.join(os.getcwd(), "output/counts/RSEM/RSEM_expected_count.tsv"),
		sample_contrast = os.path.join(os.getcwd(),config["sample_contrast"])
	run:
		# Merge RESM output files into matrix
		shell("python ./scripts/makeRSEMMatrix.py RunsByExperiment.tsv output/counts/RSEM expected_count")
		# Compute differentially expressed genes based on deg_samples.txt
		shell("run_DE_analysis.pl --matrix {params.matrix} --method DESeq2 --samples_file {params.rep_relations} --contrasts {params.sample_contrast} --output output/DEG")
		shell("cd output/counts/RSEM && PtR --matrix {params.matrix} --min_rowSums 10 -s {params.rep_relations} --log2 --CPM --sample_cor_matrix --CPM --center_rows --prin_comp 3")
		shell("cd output/DEG && analyze_diff_expr.pl --samples  {params.rep_relations} --matrix {params.matrix} -P 0.001 -C 2")
		shell("cd output/DEG && analyze_diff_expr.pl --samples {params.rep_relations} --matrix {params.matrix} -P 0.01 -C 1")

