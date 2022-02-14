#!/usr/bin/env bash

set -x
set -e

snakemake -p --rerun-incomplete --cluster-config cluster.json \
          --cluster "sbatch -N {cluster.nodes} --mem={cluster.memory} --cpus-per-task={cluster.ncpus} --parsable -A {cluster.account} -p {cluster.partition} -t {cluster.time} -o {cluster.output} -e {cluster.error}" \
          --max-jobs-per-second 50 \
          --max-status-checks-per-second 50 \
		  --jobs 50 \
		  --latency-wait 440 \
		  --notemp \
		  --cluster-status ./src/status.py \
          "$@"

