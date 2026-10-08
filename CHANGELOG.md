# Changelog

## [Unreleased]

### Added

- BBDuk pre-alignment contaminant filter (`rules/bbduk_filter.smk`). Single-pass at `k=31` removes rRNA (BBTools `ribokmers.fa.gz`), chloroplast (12 species), mitochondrion (12 species), and NCBI UniVec vector reads before STAR alignment. Toggleable via `config.yml::bbduk_enable` (off by default).
- `scripts/build_bbduk_refs.sh` — one-shot reference builder using NCBI Entrez Direct (`efetch`) plus a UniVec download from the NCBI FTP. Includes a 1-second courtesy delay between Entrez calls.
- `data/bbduk_refs/SPECIES_LIST.tsv` — canonical 12-species cp/mt accession list spanning all major land-plant clades (algae → bryophyte → lycophyte → gymnosperm → 7 angiosperm clades).
- `scripts/aggregate_bbduk_stats.py` — walks per-sample BBDuk outputs and produces a single tidy summary TSV (rRNA / cp / mt / univec / clean breakdown).
- `envs/bbduk.yaml` — conda env spec (`bioconda::bbmap=39.06`, `bioconda::entrez-direct=22.4`).
- `tests/test_bbduk_filter.sh` — smoke test that runs BBDuk on a 100k-read subset and asserts non-empty clean output + refstats.
- `strandedness` config key (`none`|`forward`|`reverse`, default `none`), passed to RSEM (`--strandedness`), featureCounts (`-s 0/1/2`) and htseq-count (`--stranded no/yes/reverse`). Previously every counter was hard-wired to unstranded. `RSEM_strandedness` is accepted as an alias.
- `scripts/create_RunsbyExperiment.py --reference <control>` writes the control second in every contrast pair. The script now also prints that a contrasts line `A<TAB>B` gives DESeq2 log2FoldChange = log2(A/B); README documents the direction under *DE contrast direction*.
- `tests/test_workflow_config.py` — static checks over the workflow sources and both config files for every fix below (none of them is reachable by `snakemake -n`).

### Changed

- `Snakefile` `align_PAIRED`, `align_SINGLE`, `alignRSEM_PAIRED`, `alignRSEM_SINGLE` now consume BBDuk clean output when `bbduk_enable: true` (backward-compatible: `false` keeps the original DAG and inputs unchanged).
- **(potentially breaking for downstream forks)** `align_SINGLE` and `alignRSEM_SINGLE` shell commands switched from positional `{input}` to named `{input.fwd_fastq}` (now that input arrives via `unpack(...)`). Forks that override these rules with custom `shell:` blocks referencing `{input}` will need to update.

### Fixed

- `alignRSEM` failed on STAR >= 2.7.11a (`EXITING: FATAL INPUT ERROR: unrecognized parameter name "quantTranscriptomeBan"`). It now uses the replacement `--quantTranscriptomeSAMoutput BanSingleEnd_BanIndels_ExtendSoftclip` (same behaviour). **STAR >= 2.7.11a is now required.** Its declared log is now actually written (`2> {log}`).
- Cluster-config resource blocks are keyed by rule name, but `config.yml` used `align_PAIRED`/`_SINGLE`, `alignRSEM_PAIRED`/`_SINGLE`, `calculateRSEMExpression_PAIRED`/`_SINGLE`, `importRaw_paired`/`_single`, `summarize`, `fastqc_raw` and `fastqc_trim_*`, none of which name a rule. STAR and RSEM jobs therefore ran on the `__default__` 2 CPU / 8 GB allocation, below `alignRSEM`'s 20 GB `--limitBAMsortRAM`. Blocks renamed to the rule names (`align`, `alignRSEM`, `calculateRSEMExpression`, `importRaw_PAIRED`, ...) in `config.yml` and `examples/example_config.yml`; keys for rules that do not exist were removed.
- `trim_PAIRED`/`trim_SINGLE` ran FastQC with `2> {log}`, overwriting fastp's report in the same log; FastQC now appends.

- `config.yml` and `examples/example_config.yml` renamed `GTF_name` → `GTFname` to match the key the Snakefile actually reads (`config["GTFname"]` / `config_dict["GTFname"]` in 6 places). Previously broke any non-RSEM read-counting run with `KeyError: 'GTFname'`.
