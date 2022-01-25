configfile: "config.json"

rule all:
    input:
        expand("bam/{sample}.bamAligned.sortedByCoord.out.bam", sample=config["samples"])

rule downloadSRA:
    input:
        expand("{reference}", reference=config["reference"])
    output:
        expand("raw_data/{sample}_{replicate}.fastq.gz",
            sample=config["samples"],
            replicate=config["replicates"])
    log:
        "logs/fastq-dump.log"
    run:
        for sample in config["samples"]:
            shell("prefetch " + sample)
            shell("parallel-fastq-dump --sra-id " + sample +
                  " --threads " + str(config["maxThreads"]) +
                  " --split-3 --outdir ./raw_data  --gzip 2> {log}")

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
