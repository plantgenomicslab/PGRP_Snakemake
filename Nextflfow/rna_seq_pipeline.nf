// RNA-seq Pipeline converted from the provided Snakemake to Nextflow (DSL2)
// Author: Friday
// Notes:
//  - Requires: STAR, samtools, fastqc, fastp, parallel-fastq-dump, rsem-calculate-expression,
//              featureCounts, htseq-count, TPMCalculator, PtR, Trinity's run_DE_analysis.pl & analyze_diff_expr.pl,
//              and the helper scripts used in the original workflow (mergeTPMCalculator.py, normalizeCounts.py,
//              calc_cdna_len.py, makeRSEMMatrix.py, summarizeNormalizedCounts.py).
//  - Configuration is supplied via nextflow.config (params.*). See the example at bottom of this file.
//  - Input sample sheet: RunsByExperiment.tsv (tab-separated) with at least columns:
//      Run, Replicate, Treatment
//    The pipeline constructs replicate-level merges and count matrices accordingly.
//  - Library layout: params.libraryType must be 'PAIRED' or 'SINGLE'.
//  - Source: params.source must be 'local' (reads in params.rawInputDir) or 'sra' (download via prefetch + parallel-fastq-dump).
//  - Read counting backends: any subset of [RSEM, featureCounts, HTseq, TPMcalculator] in params.readCounting.

nextflow.enable.dsl = 2

// ----------------------------
// Parameters & sanity checks
// ----------------------------
params.samples_file          = params.samples_file          ?: 'RunsByExperiment.tsv'
params.source                = params.source                ?: 'local'       // 'local' or 'sra'
params.libraryType           = params.libraryType           ?: 'PAIRED'      // 'PAIRED' or 'SINGLE'
params.PAIR_LIST             = params.PAIR_LIST             ?: ['_1','_2']   // used when PAIRED
params.rawInputDir           = params.rawInputDir           ?: 'raw'
params.sra_dir               = params.sra_dir               ?: 'sra'
params.genomeDir             = params.genomeDir             ?: 'genomeDir'   // STAR genomeDir path
params.ref                   = params.ref                   ?: []            // e.g., ['SA','SAindex','Genome'] to verify
params.GTFname               = params.GTFname               ?: ''
params.RSEM_prepared_genome  = params.RSEM_prepared_genome  ?: ''            // prefix path (without extension)
params.readCounting          = (params.readCounting instanceof List) ? params.readCounting : (params.readCounting ? params.readCounting.tokenize(',') : [])
params.rep_relations         = params.rep_relations         ?: 'replication_relationship.txt'
params.sample_contrast       = params.sample_contrast       ?: ''
params.runDEG                = (params.runDEG == null) ? false : params.runDEG as boolean

// Thread presets
params.threads = params.threads ?: [
  fetchSRA: 4,
  convertSRAtoFastq: 8,
  fastqc_raw: 4,
  trim: 8,
  align: 16,
  merge: 8,
  calculateRSEMExpression: 12,
  featureCount: 12,
  HTseq: 12,
  TPMCalculator: 8,
  mergeTPMCalculator: 4,
  normalizeFeatureCounts: 4,
  normalizeHTseq: 4,
  DEG_featureCounts: 8,
  DEG_HTseq: 8,
  DEG_RSEM: 8
]

// ----------------------------
// Helper: read sample sheet & build maps
// ----------------------------
process READ_SAMPLES {
  tag "parse ${params.samples_file}"
  publishDir 'output', mode: 'copy', pattern: 'replication_relationship.txt'
  input:
  val samples_file from Channel.value(params.samples_file)
  output:
  path 'parsed_samples.tsv'         // normalized rows: Run,Replicate,Treatment,Sample
  path 'replication_relationship.txt'
  val replicate_runs                into REPLICATE_RUNS
  val run_list                      into RUN_LIST
  val replicate_list                into REPLICATE_LIST

  script:
  def awk = """
    awk 'BEGIN{FS=OFS="\t"}{if(NR==1){
      for(i=1;i<=NF;i++){h[\$i]=i}
      if(!h["Run"]||!h["Replicate"]||!h["Treatment"]) {print "Missing required columns (Run, Replicate, Treatment)" > "/dev/stderr"; exit 1}
      print "Run","Replicate","Treatment","Sample"; next}
      run=\$h["Run"]; rep=\$h["Replicate"]; trt=\$h["Treatment"]; sample=run; print run,rep,trt,sample
    }' ${samples_file} > parsed_samples.tsv
  """
  def repRel = """
    awk 'BEGIN{FS=OFS="\t"}{if(NR==1){for(i=1;i<=NF;i++){h[\$i]=i}; next} print \$h["Treatment"],\$h["Replicate"]}' ${samples_file} | sort -u > replication_relationship.txt
  """
  def runsMap = """
    awk 'BEGIN{FS=OFS="\t"}{if(NR==1){for(i=1;i<=NF;i++){h[\$i]=i}; next} print \$h["Replicate"],\$h["Run"]}' ${samples_file} | sort -u \
    | awk '{a[$1]=a[$1] (a[$1]?",":"") $2} END{for(k in a) print k"\t"a[k]}' > rep_to_runs.tsv
  """
  def listRuns = "cut -f1 parsed_samples.tsv | tail -n +2 | sort -u > runs.list"
  def listReps = "cut -f2 parsed_samples.tsv | tail -n +2 | sort -u > reps.list"
  def emit = """
    echo -n "REPLICATE_RUNS\t"; cat rep_to_runs.tsv | tr '\n' ';'
  """
  def bash = """
    set -euo pipefail
    ${awk}
    ${repRel}
    ${runsMap}
    ${listRuns}
    ${listReps}
  """
  // Nextflow cannot directly emit complex maps from bash; we emit via 'val' by reading files in the nextflow layer below.
  // Thus we capture values in the when/afterScript using file reading in the executor side (Groovy block below not allowed here).
  // Workaround: we pass via file handles declared above, and at the workflow level, we re-read them.
  shell: bash
}

// Rehydrate lists/maps for channels
Channel.fromPath('reps.list').splitText().filter{ it?.trim() }.set{ REPS_CH }
Channel.fromPath('runs.list').splitText().filter{ it?.trim() }.set{ RUNS_CH }

// map of replicate -> csv of runs
def REP2RUNS = [:]
if( file('rep_to_runs.tsv').exists() ) {
  file('rep_to_runs.tsv').readLines().each { line ->
    def (rep, runsCsv) = line.split('\t')
    REP2RUNS[rep] = runsCsv.split(',') as List
  }
}

// ----------------------------
// Reference checks
// ----------------------------
process CHECK_REFERENCES {
  tag 'check refs'
  input:
  val genomeDir from Channel.value(params.genomeDir)
  val refList   from Channel.value(params.ref)
  val rsemRef   from Channel.value(params.RSEM_prepared_genome)
  output:
  val true into REF_OK
  script:
  def checks = []
  if( refList && refList.size()>0 ) {
    checks << refList.collect{ "[ -e \"${genomeDir}/${it}\" ]" }.join(' && ')
  }
  if( params.readCounting.contains('RSEM') && rsemRef ) {
    checks << "[ -e \"${rsemRef}.seq\" ]"
  }
  def cmd = checks ? checks.join(' && ') : 'true'
  """
  bash -uexo pipefail -c '${cmd}'
  """
}

// ----------------------------
// Fetch / import raw reads
// ----------------------------
process FETCH_SRA {
  tag { sample }
  when: params.source == 'sra'
  cpus params.threads.fetchSRA
  input:
  val sample from RUNS_CH
  output:
  path("${params.sra_dir}/${sample}.sra") into SRA_FILES
  script:
  """
  set -euo pipefail
  mkdir -p ${params.sra_dir}
  prefetch ${sample} --output-file ${params.sra_dir}/${sample}.sra
  vdb-validate ${params.sra_dir}/${sample}.sra
  """
}

process SRA_TO_FASTQ_PAIRED {
  tag { sample }
  when: params.source == 'sra' && params.libraryType == 'PAIRED'
  cpus params.threads.convertSRAtoFastq
  input:
  path sra from SRA_FILES.map{ it } // keep path
  output:
  tuple val(sample), path("${sample}${params.PAIR_LIST[0]}.fastq.gz"), path("${sample}${params.PAIR_LIST[1]}.fastq.gz") into RAW_PAIRED
  script:
  def sample = sra.baseName
  """
  set -euo pipefail
  parallel-fastq-dump --sra-id ${sample} --threads ${task.cpus} --split-e --gzip --outdir .
  """
}

process SRA_TO_FASTQ_SINGLE {
  tag { sample }
  when: params.source == 'sra' && params.libraryType == 'SINGLE'
  cpus params.threads.convertSRAtoFastq
  input:
  path sra from SRA_FILES
  output:
  tuple val(sample), path("${sample}.fastq.gz") into RAW_SINGLE
  script:
  def sample = sra.baseName
  """
  set -euo pipefail
  parallel-fastq-dump --sra-id ${sample} --threads ${task.cpus} --split-e --gzip --outdir .
  """
}

process IMPORT_RAW_PAIRED {
  tag { sample }
  when: params.source == 'local' && params.libraryType == 'PAIRED'
  input:
  val sample from RUNS_CH
  output:
  tuple val(sample), path("${sample}${params.PAIR_LIST[0]}.fastq.gz"), path("${sample}${params.PAIR_LIST[1]}.fastq.gz") into RAW_PAIRED
  script:
  """
  set -euo pipefail
  ln -s ${params.rawInputDir}/${sample}${params.PAIR_LIST[0]}.fastq.gz .
  ln -s ${params.rawInputDir}/${sample}${params.PAIR_LIST[1]}.fastq.gz .
  """
}

process IMPORT_RAW_SINGLE {
  tag { sample }
  when: params.source == 'local' && params.libraryType == 'SINGLE'
  input:
  val sample from RUNS_CH
  output:
  tuple val(sample), path("${sample}.fastq.gz") into RAW_SINGLE
  script:
  """
  set -euo pipefail
  ln -s ${params.rawInputDir}/${sample}.fastq.gz .
  """
}

// ----------------------------
// QC & trimming
// ----------------------------
process FASTQC_PAIRED {
  tag { sample }
  cpus params.threads.fastqc_raw
  input:
  tuple val(sample), path(r1), path(r2) from RAW_PAIRED
  output:
  tuple val(sample), path(r1), path(r2), path("${r1.baseName}_fastqc.zip"), path("${r2.baseName}_fastqc.zip") into RAW_PAIRED_QC
  script:
  """
  set -euo pipefail
  fastqc --threads ${task.cpus} ${r1} ${r2}
  """
}

process FASTQC_SINGLE {
  tag { sample }
  cpus params.threads.fastqc_raw
  input:
  tuple val(sample), path(r) from RAW_SINGLE
  output:
  tuple val(sample), path(r), path("${r.baseName}_fastqc.zip") into RAW_SINGLE_QC
  script:
  """
  set -euo pipefail
  fastqc --threads ${task.cpus} ${r}
  """
}

process TRIM_PAIRED {
  tag { sample }
  cpus params.threads.trim
  input:
  tuple val(sample), path(r1), path(r2), path(q1), path(q2) from RAW_PAIRED_QC
  output:
  tuple val(sample), path("${sample}${params.PAIR_LIST[0]}_trimmed.fq.gz"), path("${sample}${params.PAIR_LIST[1]}_trimmed.fq.gz") into TRIMMED_PAIRED
  script:
  """
  set -euo pipefail
  fastp --detect_adapter_for_pe --overrepresentation_analysis --cut_right --thread ${task.cpus} \
    --html ${sample}.fastp.html --json ${sample}.fastp.json \
    -i ${r1} -I ${r2} -o ${sample}${params.PAIR_LIST[0]}_trimmed.fq.gz -O ${sample}${params.PAIR_LIST[1]}_trimmed.fq.gz
  fastqc --threads ${task.cpus} ${sample}${params.PAIR_LIST[0]}_trimmed.fq.gz ${sample}${params.PAIR_LIST[1]}_trimmed.fq.gz
  """
}

process TRIM_SINGLE {
  tag { sample }
  cpus params.threads.trim
  input:
  tuple val(sample), path(r), path(qc) from RAW_SINGLE_QC
  output:
  tuple val(sample), path("${sample}_trimmed.fq.gz") into TRIMMED_SINGLE
  script:
  """
  set -euo pipefail
  fastp --overrepresentation_analysis --cut_right --thread ${task.cpus} \
    --html ${sample}.fastp.html --json ${sample}.fastp.json \
    -i ${r} -o ${sample}_trimmed.fq.gz
  fastqc --threads ${task.cpus} ${sample}_trimmed.fq.gz
  """
}

// ----------------------------
// STAR alignment (genome BAM) & STAR-to-transcriptome (RSEM)
// ----------------------------
process STAR_ALIGN_PAIRED {
  tag { sample }
  cpus params.threads.align
  input:
  tuple val(sample), path(t1), path(t2) from TRIMMED_PAIRED
  output:
  tuple val(sample), path("${sample}.bamAligned.sortedByCoord.out.bam") into STAR_BAMS
  script:
  """
  set -euo pipefail
  STAR --runMode alignReads --runThreadN ${task.cpus} \
    --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 \
    --quantMode TranscriptomeSAM GeneCounts \
    --outBAMsortingBinsN 200 \
    --genomeDir ${params.genomeDir} \
    --readFilesCommand gunzip -c --readFilesIn ${t1} ${t2} \
    --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}.bam
  samtools flagstat ${sample}.bamAligned.sortedByCoord.out.bam
  """
}

process STAR_ALIGN_SINGLE {
  tag { sample }
  cpus params.threads.align
  input:
  tuple val(sample), path(t) from TRIMMED_SINGLE
  output:
  tuple val(sample), path("${sample}.bamAligned.sortedByCoord.out.bam") into STAR_BAMS
  script:
  """
  set -euo pipefail
  STAR --runMode alignReads --runThreadN ${task.cpus} \
    --outFilterMultimapNmax 100 --alignIntronMin 25 --alignIntronMax 50000 \
    --quantMode TranscriptomeSAM GeneCounts \
    --outBAMsortingBinsN 200 \
    --genomeDir ${params.genomeDir} \
    --readFilesCommand gunzip -c --readFilesIn ${t} \
    --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample}.bam
  samtools flagstat ${sample}.bamAligned.sortedByCoord.out.bam
  """
}

process STAR_ALIGN_RSEM_PAIRED {
  tag { sample }
  when: params.readCounting.contains('RSEM')
  cpus params.threads.align
  input:
  tuple val(sample), path(t1), path(t2) from TRIMMED_PAIRED
  output:
  tuple val(sample), path("${sample}.xs.bamAligned.toTranscriptome.out.bam") into RSEM_BAMS
  script:
  """
  set -euo pipefail
  STAR --runMode alignReads --runThreadN ${task.cpus} \
    --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 25000 \
    --limitBAMsortRAM 20000000000 --outBAMsortingBinsN 200 \
    --genomeDir ${params.genomeDir} \
    --readFilesCommand gunzip -c --readFilesIn ${t1} ${t2} \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
    --quantTranscriptomeBan IndelSoftclipSingleend --alignEndsType EndToEnd \
    --outFileNamePrefix ${sample}.xs.bam
  """
}

process STAR_ALIGN_RSEM_SINGLE {
  tag { sample }
  when: params.readCounting.contains('RSEM')
  cpus params.threads.align
  input:
  tuple val(sample), path(t) from TRIMMED_SINGLE
  output:
  tuple val(sample), path("${sample}.xs.bamAligned.toTranscriptome.out.bam") into RSEM_BAMS
  script:
  """
  set -euo pipefail
  STAR --runMode alignReads --runThreadN ${task.cpus} \
    --outFilterMultimapNmax 10 --alignIntronMin 25 --alignIntronMax 25000 \
    --limitBAMsortRAM 20000000000 --outBAMsortingBinsN 200 \
    --genomeDir ${params.genomeDir} \
    --readFilesCommand gunzip -c --readFilesIn ${t} \
    --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM \
    --quantTranscriptomeBan IndelSoftclipSingleend --alignEndsType EndToEnd \
    --outFileNamePrefix ${sample}.xs.bam
  """
}

// ----------------------------
// Merge runs per replicate (STAR genome BAMs and RSEM transcriptome BAMs)
// ----------------------------
process MERGE_STAR_BY_REP {
  tag { rep }
  cpus params.threads.merge
  input:
  val rep from REPS_CH
  tuple val(sample), path(bam) from STAR_BAMS.collect().flatten().filter{ it instanceof Tuple }.map{ it }
  // We'll select the BAMs for this rep below
  output:
  tuple val(rep), path("${rep}.STAR.bam") into STAR_BY_REP
  script:
  def runs = REP2RUNS[rep] ?: []
  def inList = runs.collect{ "${it}.bamAligned.sortedByCoord.out.bam" }.findAll{ file(it).exists() }
  if( inList.size() == 0 ) {
    return "echo No BAMs for ${rep} > ${rep}.STAR.bam && false" // fail fast if mismatch
  }
  def mergeCmd = (inList.size()>1) ? "samtools merge -@ ${task.cpus} ${rep}.STAR.bam ${inList.join(' ')} && samtools index ${rep}.STAR.bam" : "ln -s ${inList[0]} ${rep}.STAR.bam"
  """
  set -euo pipefail
  ${mergeCmd}
  samtools flagstat ${rep}.STAR.bam || true
  """
}

process MERGE_RSEM_BY_REP {
  tag { rep }
  when: params.readCounting.contains('RSEM')
  cpus params.threads.merge
  input:
  val rep from REPS_CH
  output:
  tuple val(rep), path("${rep}.RSEM.bam") into RSEM_BY_REP
  script:
  def runs = REP2RUNS[rep] ?: []
  def inList = runs.collect{ "${it}.xs.bamAligned.toTranscriptome.out.bam" }.findAll{ file(it).exists() }
  if( inList.size() == 0 ) {
    return "echo No RSEM BAMs for ${rep} > ${rep}.RSEM.bam && false"
  }
  def mergeCmd = (inList.size()>1) ? "samtools merge -@ ${task.cpus} ${rep}.RSEM.bam ${inList.join(' ')} && samtools index ${rep}.RSEM.bam" : "ln -s ${inList[0]} ${rep}.RSEM.bam"
  """
  set -euo pipefail
  ${mergeCmd}
  samtools flagstat ${rep}.RSEM.bam || true
  """
}

// ----------------------------
// Counting backends
// ----------------------------
process RSEM_CALC {
  tag { rep }
  when: params.readCounting.contains('RSEM')
  cpus params.threads.calculateRSEMExpression
  input:
  tuple val(rep), path(bam) from RSEM_BY_REP
  output:
  path("${rep}.genes.results") into RSEM_RESULTS
  script:
  """
  set -euo pipefail
  rsem-calculate-expression --bam --no-bam-output -p ${task.cpus} ${bam} ${params.RSEM_prepared_genome} ${rep}
  mv ${rep}.genes.results ${rep}.genes.results
  """
}

process FEATURECOUNTS {
  tag 'featureCounts'
  when: params.readCounting.contains('featureCounts')
  cpus params.threads.featureCount
  input:
  tuple val(rep), path(bam) from STAR_BY_REP
  // need all bams together; use collectFile pattern
  // We build a file list and run once
  output:
  path 'featureCounts.cnt' into FC_RAW
  path 'featureCount_clean.cnt' into FC_CLEAN
  script:
  def bams = Channel.fromPath('*.STAR.bam').getVal() // not allowed in script; alternative: run featureCounts on glob
  """
  set -euo pipefail
  featureCounts -o featureCounts.cnt -T ${task.cpus} -Q 1 -p -M -g gene_id -a ${params.GTFname} *.STAR.bam 2> featureCounts.log
  cat featureCounts.cnt | egrep -v '#' | sed 's/\\Aligned\.sortedByCoord\.out\.bam//g; s/\.bam//g; s/.*\///g' > featureCount_clean.cnt
  """
}

process HTSEQ_COUNT {
  tag 'HTSeq'
  when: params.readCounting.contains('HTseq')
  cpus params.threads.HTseq
  input:
  tuple val(rep), path(bam) from STAR_BY_REP
  output:
  path 'htseq-count.tsv' into HTSEQ_TSV
  script:
  """
  set -euo pipefail
  htseq-count --format bam --order pos --mode union --stranded=no --type exon --idattr gene_id --nprocesses ${task.cpus} --counts_output htseq-count.tsv *.STAR.bam ${params.GTFname} &> HTseq.log
  sed -i '1 i\gene\t'htseq-count.tsv
  sed -i 's/\.bam//g' htseq-count.tsv
  """
}

process TPMCALC {
  tag { rep }
  when: params.readCounting.contains('TPMcalculator')
  cpus params.threads.TPMCalculator
  input:
  tuple val(rep), path(bam) from STAR_BY_REP
  output:
  path("${rep}_genes.out") into TPM_OUT_PER_REP
  script:
  """
  set -euo pipefail
  TPMCalculator -k gene_id -t transcript_id -o 0 -g ${params.GTFname} -b ${bam}
  mv ${rep}_genes.* .
  """
}

process MERGE_TPMCALC {
  tag 'merge TPMcalc'
  when: params.readCounting.contains('TPMcalculator')
  cpus params.threads.mergeTPMCalculator
  input:
  path outs from TPM_OUT_PER_REP.collect()
  output:
  path 'tpmcalculator-merged.tsv' into TPM_MERGED
  script:
  """
  set -euo pipefail
  ./scripts/mergeTPMCalculator.py ${outs.join(' ')}
  """
}

// ----------------------------
// Normalization & DEG (optional)
// ----------------------------
process NORMALIZE_FEATURECOUNTS {
  when: params.readCounting.contains('featureCounts')
  cpus params.threads.normalizeFeatureCounts
  input:
  path cnt from FC_CLEAN
  output:
  path 'featureCounts.tpm.tsv' into FC_TPM
  path 'featureCounts.fpkm.tsv' into FC_FPKM
  script:
  """
  set -euo pipefail
  python ./scripts/normalizeCounts.py featureCounts ${params.GTFname} ${cnt} featureCounts ${params.genomeDir}
  """
}

process NORMALIZE_HTSEQ {
  when: params.readCounting.contains('HTseq')
  cpus params.threads.normalizeHTseq
  input:
  path tsv from HTSEQ_TSV
  output:
  path 'htseq-count.tpm.tsv' into HTSEQ_TPM
  path 'htseq-count.fpkm.tsv' into HTSEQ_FPKM
  script:
  """
  set -euo pipefail
  ./scripts/calc_cdna_len.py ${params.GTFname} gene_id > ${params.genomeDir}/cds_length.tsv
  python scripts/normalizeCounts.py HTseq ${params.GTFname} ${tsv} htseq-count ${params.genomeDir}
  """
}

process SUMMARIZE_TPMCALC {
  when: params.readCounting.contains('TPMcalculator')
  input:
  path m from TPM_MERGED
  output:
  path 'tpmcalculator-merged.tsv.average.tsv'
  script:
  """
  ./scripts/summarizeNormalizedCounts.py Gene_Id ${m}
  """
}

process SUMMARIZE_HTSEQ {
  when: params.readCounting.contains('HTseq')
  input:
  path tpm from HTSEQ_TPM
  path fpkm from HTSEQ_FPKM
  output:
  path 'htseq-count.tpm.tsv.average.tsv'
  script:
  """
  ./scripts/summarizeNormalizedCounts.py gene ${tpm}
  ./scripts/summarizeNormalizedCounts.py gene ${fpkm}
  """
}

process SUMMARIZE_FEATURECOUNTS {
  when: params.readCounting.contains('featureCounts')
  input:
  path tpm from FC_TPM
  path fpkm from FC_FPKM
  output:
  path 'featureCounts.tpm.tsv.average.tsv'
  script:
  """
  ./scripts/summarizeNormalizedCounts.py Geneid ${tpm}
  ./scripts/summarizeNormalizedCounts.py Geneid ${fpkm}
  """
}

process SUMMARIZE_RSEM {
  when: params.readCounting.contains('RSEM')
  input:
  path results from RSEM_RESULTS.collect()
  output:
  path 'RSEM_TPM.tsv'
  path 'RSEM_TPM.tsv.average.tsv'
  script:
  """
  python ./scripts/makeRSEMMatrix.py ${params.samples_file} . TPM
  python ./scripts/makeRSEMMatrix.py ${params.samples_file} . FPKM
  ./scripts/summarizeNormalizedCounts.py gene_id RSEM_TPM.tsv
  ./scripts/summarizeNormalizedCounts.py gene_id RSEM_FPKM.tsv
  """
}

// DEG: featureCounts
process DEG_FEATURECOUNTS {
  when: params.runDEG && params.readCounting.contains('featureCounts') && params.sample_contrast
  cpus params.threads.DEG_featureCounts
  input:
  path cnt from FC_CLEAN
  output:
  path "featureCount_clean.cnt.${params.sample_contrast.tokenize('\n')? 'DESeq2.DE_results' : 'DESeq2.DE_results'}" optional true
  script:
  """
  set -euo pipefail
  run_DE_analysis.pl --matrix ${cnt} --method DESeq2 --samples_file ${params.rep_relations} --contrasts ${params.sample_contrast} --output .
  analyze_diff_expr.pl --samples ${params.rep_relations} --matrix ${cnt} -P 0.001 -C 2
  analyze_diff_expr.pl --samples ${params.rep_relations} --matrix ${cnt} -P 0.01 -C 1
  PtR --matrix featureCount_clean.cnt --min_rowSums 10 -s ${params.rep_relations} --log2 --CPM --sample_cor_matrix --CPM --center_rows --prin_comp 3
  """
}

// DEG: HTSeq
process DEG_HTSEQ {
  when: params.runDEG && params.readCounting.contains('HTseq') && params.sample_contrast
  cpus params.threads.DEG_HTseq
  input:
  path tsv from HTSEQ_TSV
  script:
  """
  set -euo pipefail
  run_DE_analysis.pl --matrix ${tsv} --method DESeq2 --samples_file ${params.rep_relations} --contrasts ${params.sample_contrast} --output .
  analyze_diff_expr.pl --samples ${params.rep_relations} --matrix ${tsv} -P 0.001 -C 2
  analyze_diff_expr.pl --samples ${params.rep_relations} --matrix ${tsv} -P 0.01 -C 1
  PtR --matrix htseq-count.tsv --min_rowSums 10 -s ${params.rep_relations} --log2 --CPM --sample_cor_matrix --CPM --center_rows --prin_comp 3
  """
}

// DEG: RSEM
process DEG_RSEM {
  when: params.runDEG && params.readCounting.contains('RSEM') && params.sample_contrast
  cpus params.threads.DEG_RSEM
  input:
  path rs from RSEM_RESULTS.collect()
  script:
  """
  set -euo pipefail
  python ./scripts/makeRSEMMatrix.py ${params.samples_file} . expected_count
  run_DE_analysis.pl --matrix RSEM_expected_count.tsv --method DESeq2 --samples_file ${params.rep_relations} --contrasts ${params.sample_contrast} --output .
  PtR --matrix RSEM_expected_count.tsv --min_rowSums 10 -s ${params.rep_relations} --log2 --CPM --sample_cor_matrix --CPM --center_rows --prin_comp 3
  analyze_diff_expr.pl --samples ${params.rep_relations} --matrix RSEM_expected_count.tsv -P 0.001 -C 2
  analyze_diff_expr.pl --samples ${params.rep_relations} --matrix RSEM_expected_count.tsv -P 0.01 -C 1
  """
}

// ----------------------------
// Workflow definition
// ----------------------------
workflow {
  REF_OK
  READ_SAMPLES()

  if( params.libraryType == 'PAIRED' ) {
    if( params.source == 'sra' ) SRA_TO_FASTQ_PAIRED(FETCH_SRA)
    else                         IMPORT_RAW_PAIRED()
    FASTQC_PAIRED(TRIM_PAIRED(RAW_PAIRED_QC(FASTQC_PAIRED))) // chaining ensures order
    STAR_ALIGN_PAIRED(TRIM_PAIRED)
    if( params.readCounting.contains('RSEM') ) STAR_ALIGN_RSEM_PAIRED(TRIM_PAIRED)
  } else {
    if( params.source == 'sra' ) SRA_TO_FASTQ_SINGLE(FETCH_SRA)
    else                         IMPORT_RAW_SINGLE()
    FASTQC_SINGLE(TRIM_SINGLE(RAW_SINGLE_QC(FASTQC_SINGLE)))
    STAR_ALIGN_SINGLE(TRIM_SINGLE)
    if( params.readCounting.contains('RSEM') ) STAR_ALIGN_RSEM_SINGLE(TRIM_SINGLE)
  }

  MERGE_STAR_BY_REP()
  if( params.readCounting.contains('RSEM') ) MERGE_RSEM_BY_REP()

  if( params.readCounting.contains('RSEM') ) {
    RSEM_CALC()
    SUMMARIZE_RSEM(RSEM_RESULTS)
    if( params.runDEG && params.sample_contrast ) DEG_RSEM(RSEM_RESULTS)
  }

  if( params.readCounting.contains('featureCounts') ) {
    FEATURECOUNTS(STAR_BY_REP)
    NORMALIZE_FEATURECOUNTS(FC_CLEAN)
    SUMMARIZE_FEATURECOUNTS(FC_TPM, FC_FPKM)
    if( params.runDEG && params.sample_contrast ) DEG_FEATURECOUNTS(FC_CLEAN)
  }

  if( params.readCounting.contains('HTseq') ) {
    HTSEQ_COUNT(STAR_BY_REP)
    NORMALIZE_HTSEQ(HTSEQ_TSV)
    SUMMARIZE_HTSEQ(HTSEQ_TPM, HTSEQ_FPKM)
    if( params.runDEG && params.sample_contrast ) DEG_HTSEQ(HTSEQ_TSV)
  }

  if( params.readCounting.contains('TPMcalculator') ) {
    TPMCALC(STAR_BY_REP)
    MERGE_TPMCALC(TPM_OUT_PER_REP)
    SUMMARIZE_TPMCALC(TPM_MERGED)
  }
}

// ----------------------------
// Example nextflow.config (place in nextflow.config)
// ----------------------------
/*
params {
  samples_file = 'RunsByExperiment.tsv'
  source       = 'local'   // or 'sra'
  libraryType  = 'PAIRED'  // or 'SINGLE'
  PAIR_LIST    = ['_1','_2']
  rawInputDir  = '/path/to/fastq'
  sra_dir      = 'sra'
  genomeDir    = '/path/to/STAR_genomeDir/'
  ref          = ['Genome','SA','SAindex']
  GTFname      = '/path/to/genes.gtf'
  RSEM_prepared_genome = '/path/to/rsem/ref_prefix'
  readCounting = ['RSEM','featureCounts','HTseq','TPMcalculator']
  rep_relations = 'replication_relationship.txt'
  sample_contrast = 'sample_contrasts.txt' // two-column tsv, e.g., Top\tBottom
  runDEG = true

  threads = [
    fetchSRA: 4, convertSRAtoFastq: 8, fastqc_raw: 4, trim: 8, align: 16, merge: 8,
    calculateRSEMExpression: 12, featureCount: 12, HTseq: 12, TPMCalculator: 8, mergeTPMCalculator: 4,
    normalizeFeatureCounts: 4, normalizeHTseq: 4, DEG_featureCounts: 8, DEG_HTseq: 8, DEG_RSEM: 8
  ]
}

process.executor = 'local'
process.withName:STAR_ALIGN_* { cpus = params.threads.align; memory = '32 GB' }
process.withName:TRIM_*       { cpus = params.threads.trim;  memory = '8 GB' }
process.withName:FEATURECOUNTS{ cpus = params.threads.featureCount }
// Add container or conda profiles as needed
*/


// =============================================================
// Separate files: params.yaml and resources-slurm.config
// =============================================================

// ===== File: params.yaml =====
// Dataset- and behavior-level parameters only (no resources here)
// Load with:  -params-file params.yaml
# Samples
samples_file: RunsByExperiment.tsv
source: local            # local | sra
libraryType: PAIRED      # PAIRED | SINGLE
PAIR_LIST: ['_1','_2']
rawInputDir: /path/to/fastq
sra_dir: sra

# References
genomeDir: /path/to/STAR_genomeDir
ref: ['Genome','SA','SAindex']      # files expected under genomeDir
GTFname: /path/to/genes.gtf
RSEM_prepared_genome: /path/to/rsem/ref_prefix

# Counting backends
readCounting: ['RSEM','featureCounts','HTseq','TPMcalculator']

# DEG options
rep_relations: replication_relationship.txt
sample_contrast: sample_contrasts.txt  # two-column TSV (Top	Bottom)
runDEG: true


// ===== File: resources-slurm.config =====
// Resource & executor presets for Slurm clusters
// Load with:  -c resources-slurm.config
process {
  executor       = 'slurm'
  queue          = 'batch'          // change to your partition
  errorStrategy  = 'retry'
  maxRetries     = 2
  maxForks       = 50               // throttle concurrent jobs

  // Per-process resource overrides (these override cpus/mem in main.nf)
  withName: READ_SAMPLES          { cpus = 1;  memory = '1 GB';  time = '10m' }
  withName: CHECK_REFERENCES      { cpus = 1;  memory = '1 GB';  time = '10m' }

  withName: FETCH_SRA             { cpus = 4;  memory = '8 GB';  time = '2h' }
  withName: SRA_TO_FASTQ_PAIRED   { cpus = 8;  memory = '16 GB'; time = '4h' }
  withName: SRA_TO_FASTQ_SINGLE   { cpus = 8;  memory = '16 GB'; time = '4h' }
  withName: IMPORT_RAW_*          { cpus = 1;  memory = '1 GB';  time = '30m' }

  withName: FASTQC_*              { cpus = 4;  memory = '8 GB';  time = '1h' }
  withName: TRIM_*                { cpus = 8;  memory = '16 GB'; time = '2h' }

  withName: STAR_ALIGN_*          { cpus = 16; memory = '32 GB'; time = '12h' }
  withName: MERGE_STAR_BY_REP     { cpus = 8;  memory = '16 GB'; time = '2h' }

  withName: STAR_ALIGN_RSEM_*     { cpus = 16; memory = '32 GB'; time = '12h' }
  withName: MERGE_RSEM_BY_REP     { cpus = 8;  memory = '16 GB'; time = '2h' }
  withName: RSEM_CALC             { cpus = 12; memory = '24 GB'; time = '8h' }

  withName: FEATURECOUNTS         { cpus = 12; memory = '16 GB'; time = '4h' }
  withName: HTSEQ_COUNT           { cpus = 12; memory = '24 GB'; time = '6h' }
  withName: TPMCALC               { cpus = 8;  memory = '16 GB'; time = '3h' }
  withName: MERGE_TPMCALC         { cpus = 4;  memory = '8 GB';  time = '1h' }

  withName: NORMALIZE_*           { cpus = 4;  memory = '8 GB';  time = '1h' }
  withName: SUMMARIZE_*           { cpus = 2;  memory = '4 GB';  time = '30m' }

  withName: DEG_*                 { cpus = 8;  memory = '32 GB'; time = '8h' }

  // Optional: add Slurm extras (uncomment & edit as needed)
  // clusterOptions = '--account=YOUR_ACC --qos=YOUR_QOS'
}

executor {
  name = 'slurm'
  submitRateLimit = '10 sec'
}

// (Optional) built-in reports can also be enabled via CLI flags instead of config
// timeline { enabled = true }
// report   { enabled = true }
// trace    { enabled = true }
// dag      { enabled = true }


// ===== File: nf.sbatch (optional wrapper to submit Nextflow itself to Slurm) =====
// Usage:  sbatch nf.sbatch
#!/usr/bin/env bash
#SBATCH -J rnaseq_nf
#SBATCH -p batch                 # change to your partition
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 48:00:00
#SBATCH -o nf.%j.out
#SBATCH -e nf.%j.err

module purge
# module load nextflow/24.10  # or your site module
# module load singularity     # if using containers

# JVM options for the Nextflow launcher (not for tasks)
export NXF_OPTS="-Xms1g -Xmx4g"

nextflow run main.nf \
  -params-file params.yaml \
  -c resources-slurm.config \
  -resume \
  -with-report report.html \
  -with-trace trace.txt \
  -with-timeline timeline.html \
  -with-dag flowchart.png


// ===== Quick commands =====
// From a login node (recommended on many clusters):
//   nextflow run main.nf -params-file params.yaml -c resources-slurm.config -resume \
//     -with-report report.html -with-trace trace.txt -with-timeline timeline.html -with-dag flowchart.png
// Or via the wrapper:
//   sbatch nf.sbatch
