# Create the two files requested and make them available for download
params_yaml = """# params.yaml — dataset & behavior settings (no resources here)
# Load with: -params-file params.yaml

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
sample_contrast: sample_contrasts.txt   # two-column TSV (Top\tBottom)
runDEG: true
"""

resources_slurm = """// resources-slurm.config — process resources & Slurm executor
// Load with: -c resources-slurm.config

process {
  executor       = 'slurm'
  queue          = 'batch'          // change to your partition
  errorStrategy  = 'retry'
  maxRetries     = 2
  maxForks       = 50               // throttle concurrent jobs

  // Lightweight helpers
  withName: 'READ_SAMPLES'      { cpus = 1;  memory = '1 GB';  time = '10m' }
  withName: 'CHECK_REFERENCES'  { cpus = 1;  memory = '1 GB';  time = '10m' }

  // Ingest
  withName: 'FETCH_SRA'         { cpus = 4;  memory = '8 GB';  time = '2h' }
  withName: /SRA_TO_FASTQ_.*/   { cpus = 8;  memory = '16 GB'; time = '4h' }
  withName: /IMPORT_RAW_.*/     { cpus = 1;  memory = '1 GB';  time = '30m' }

  // QC & trimming
  withName: /FASTQC_.*/         { cpus = 4;  memory = '8 GB';  time = '1h' }
  withName: /TRIM_.*/           { cpus = 8;  memory = '16 GB'; time = '2h' }

  // Alignment & merges
  withName: /STAR_ALIGN_.*/     { cpus = 16; memory = '32 GB'; time = '12h' }
  withName: 'MERGE_STAR_BY_REP' { cpus = 8;  memory = '16 GB'; time = '2h' }

  // RSEM path
  withName: /STAR_ALIGN_RSEM_.*/{ cpus = 16; memory = '32 GB'; time = '12h' }
  withName: 'MERGE_RSEM_BY_REP' { cpus = 8;  memory = '16 GB'; time = '2h' }
  withName: 'RSEM_CALC'         { cpus = 12; memory = '24 GB'; time = '8h' }

  // Counting backends
  withName: 'FEATURECOUNTS'     { cpus = 12; memory = '16 GB'; time = '4h' }
  withName: 'HTSEQ_COUNT'       { cpus = 12; memory = '24 GB'; time = '6h' }
  withName: 'TPMCALC'           { cpus = 8;  memory = '16 GB'; time = '3h' }
  withName: 'MERGE_TPMCALC'     { cpus = 4;  memory = '8 GB';  time = '1h' }

  // Normalization/Summary
  withName: /NORMALIZE_.*/      { cpus = 4;  memory = '8 GB';  time = '1h' }
  withName: /SUMMARIZE_.*/      { cpus = 2;  memory = '4 GB';  time = '30m' }

  // DEG
  withName: /DEG_.*/            { cpus = 8;  memory = '32 GB'; time = '8h' }

  // Optional Slurm flags (uncomment & edit)
  // clusterOptions = '--account=YOUR_ACC --qos=YOUR_QOS'
}

executor {
  name = 'slurm'
  submitRateLimit = '10 sec'
}
"""

with open('/mnt/data/params.yaml', 'w') as f:
    f.write(params_yaml)

with open('/mnt/data/resources-slurm.config', 'w') as f:
    f.write(resources_slurm)

print("Files created:\n- /mnt/data/params.yaml\n- /mnt/data/resources-slurm.config")
