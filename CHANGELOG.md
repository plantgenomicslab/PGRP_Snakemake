# Changelog

## [Unreleased]

### Added

- BBDuk pre-alignment contaminant filter (`rules/bbduk_filter.smk`). Single-pass at `k=31` removes rRNA (BBTools `ribokmers.fa.gz`), chloroplast (12 species), mitochondrion (12 species), and NCBI UniVec vector reads before STAR alignment. Toggleable via `config.yml::bbduk_enable` (off by default).
- `scripts/build_bbduk_refs.sh` — one-shot reference builder using NCBI Entrez Direct (`efetch`) plus a UniVec download from the NCBI FTP. Includes a 1-second courtesy delay between Entrez calls.
- `data/bbduk_refs/SPECIES_LIST.tsv` — canonical 12-species cp/mt accession list spanning all major land-plant clades (algae → bryophyte → lycophyte → gymnosperm → 7 angiosperm clades).
- `scripts/aggregate_bbduk_stats.py` — walks per-sample BBDuk outputs and produces a single tidy summary TSV (rRNA / cp / mt / univec / clean breakdown).
- `envs/bbduk.yaml` — conda env spec (`bioconda::bbmap=39.06`, `bioconda::entrez-direct=22.4`).
- `tests/test_bbduk_filter.sh` — smoke test that runs BBDuk on a 100k-read subset and asserts non-empty clean output + refstats.

### Changed

- `Snakefile` `align_PAIRED`, `align_SINGLE`, `alignRSEM_PAIRED`, `alignRSEM_SINGLE` now consume BBDuk clean output when `bbduk_enable: true` (backward-compatible: `false` keeps the original DAG and inputs unchanged).
- **(potentially breaking for downstream forks)** `align_SINGLE` and `alignRSEM_SINGLE` shell commands switched from positional `{input}` to named `{input.fwd_fastq}` (now that input arrives via `unpack(...)`). Forks that override these rules with custom `shell:` blocks referencing `{input}` will need to update.

### Fixed

- `config.yml` and `examples/example_config.yml` renamed `GTF_name` → `GTFname` to match the key the Snakefile actually reads (`config["GTFname"]` / `config_dict["GTFname"]` in 6 places). Previously broke any non-RSEM read-counting run with `KeyError: 'GTFname'`.
