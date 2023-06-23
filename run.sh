#!/usr/bin/env bash

set -x
set -e

snakemake -p --rerun-incomplete --cluster-config config.yml \
		--snakefile Snakefile \
                --max-jobs-per-second 50 \
		--max-status-checks-per-second 50 \
		--jobs 150 \
		--latency-wait 30 \
		--notemp \
		--cluster-status ./scripts/status.py \
		--cluster "sbatch -N {cluster.nodes} --mem={cluster.memory} --cpus-per-task={cluster.ncpus} \
				-J {cluster.name} \
				--parsable -A {cluster.account} -p {cluster.partition} \
				-t {cluster.time} -o {cluster.output} -e {cluster.error}" \
		"$@"

