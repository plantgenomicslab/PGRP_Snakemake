# BBDuk pre-alignment contaminant filter.
# References: rRNA k-mers (BBTools ribokmers.fa.gz) + 12-species cp + 12-species mt + UniVec.
# Single-pass at k=31; outputs clean (-> STAR) + contam (archived) + per-sample stats.

import os
import shutil
import glob


def _resolve_ribokmers():
    """Locate ribokmers.fa.gz inside the BBTools install.

    'auto' searches relative to $(which bbduk.sh); otherwise the config value
    is treated as an explicit path.
    """
    cfg = config.get("bbduk_ribokmers", "auto")
    if cfg and cfg != "auto":
        if not os.path.isfile(cfg):
            raise FileNotFoundError(
                f"bbduk_ribokmers={cfg} does not exist. Set 'auto' or a valid path."
            )
        return cfg
    bbduk = shutil.which(config.get("bbduk_bin", "bbduk.sh"))
    if not bbduk:
        raise FileNotFoundError(
            "bbduk.sh not on PATH; activate the bbmap conda env (envs/bbduk.yaml) "
            "or set bbduk_bin to an absolute path in config.yml."
        )
    base = os.path.dirname(bbduk)
    candidates = []
    for layout in ("opt/bbmap-*", "share/bbmap-*"):
        candidates += glob.glob(os.path.join(base, "..", layout, "resources", "ribokmers.fa.gz"))
    if not candidates:
        raise FileNotFoundError(
            "ribokmers.fa.gz not found near bbduk.sh; set bbduk_ribokmers explicitly in config.yml."
        )
    return candidates[0]


def _bbduk_ref_string(*_args, **_kwargs):
    """Lazy ref-string builder. Resolved at job-execution time, not parse time,
    so bbduk.sh need only be on PATH when the rule actually runs."""
    refs_dir = config["bbduk_refs_dir"]
    parts = [
        _resolve_ribokmers(),
        os.path.join(refs_dir, "univec.fa.gz"),
        os.path.join(refs_dir, "cp_12sp.fa.gz"),
        os.path.join(refs_dir, "mt_12sp.fa.gz"),
    ]
    return ",".join(parts)


rule bbduk_filter_PAIRED:
    input:
        fwd_fastq = "output/{replicate}/{sample}/trim/{sample}" + PAIR_LIST[0] + "_trimmed.fq.gz",
        rev_fastq = "output/{replicate}/{sample}/trim/{sample}" + PAIR_LIST[1] + "_trimmed.fq.gz",
    output:
        fwd_clean = "output/{replicate}/{sample}/bbduk/{sample}" + PAIR_LIST[0] + "_clean.fq.gz",
        rev_clean = "output/{replicate}/{sample}/bbduk/{sample}" + PAIR_LIST[1] + "_clean.fq.gz",
        fwd_contam = "output/{replicate}/{sample}/bbduk/{sample}" + PAIR_LIST[0] + "_contam.fq.gz",
        rev_contam = "output/{replicate}/{sample}/bbduk/{sample}" + PAIR_LIST[1] + "_contam.fq.gz",
        stats    = "output/{replicate}/{sample}/bbduk/{sample}_stats.txt",
        refstats = "output/{replicate}/{sample}/bbduk/{sample}_refstats.txt",
    log: "output/{replicate}/{sample}/logs/{sample}_bbduk.log"
    threads: config["bbduk_threads"]
    resources:
        mem_mb = config["bbduk_mem_gb"] * 1024,
    params:
        bin   = config["bbduk_bin"],
        refs  = _bbduk_ref_string,
        k     = config["bbduk_k"],
        mcf   = config["bbduk_mcf"],
        mem_g = config["bbduk_mem_gb"],
        extra = config.get("bbduk_extra_args", ""),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fwd_clean})"
        {params.bin} \
            in1={input.fwd_fastq} in2={input.rev_fastq} \
            out1={output.fwd_clean} out2={output.rev_clean} \
            outm1={output.fwd_contam} outm2={output.rev_contam} \
            ref={params.refs} \
            k={params.k} mcf={params.mcf} \
            stats={output.stats} refstats={output.refstats} \
            threads={threads} -Xmx{params.mem_g}g \
            {params.extra} 2> {log}
        """


rule bbduk_filter_SINGLE:
    input:
        fwd_fastq = "output/{replicate}/{sample}/trim/{sample}_trimmed.fq.gz",
    output:
        fwd_clean  = "output/{replicate}/{sample}/bbduk/{sample}_clean.fq.gz",
        fwd_contam = "output/{replicate}/{sample}/bbduk/{sample}_contam.fq.gz",
        stats    = "output/{replicate}/{sample}/bbduk/{sample}_stats.txt",
        refstats = "output/{replicate}/{sample}/bbduk/{sample}_refstats.txt",
    log: "output/{replicate}/{sample}/logs/{sample}_bbduk.log"
    threads: config["bbduk_threads"]
    resources:
        mem_mb = config["bbduk_mem_gb"] * 1024,
    params:
        bin   = config["bbduk_bin"],
        refs  = _bbduk_ref_string,
        k     = config["bbduk_k"],
        mcf   = config["bbduk_mcf"],
        mem_g = config["bbduk_mem_gb"],
        extra = config.get("bbduk_extra_args", ""),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fwd_clean})"
        {params.bin} \
            in={input.fwd_fastq} \
            out={output.fwd_clean} \
            outm={output.fwd_contam} \
            ref={params.refs} \
            k={params.k} mcf={params.mcf} \
            stats={output.stats} refstats={output.refstats} \
            threads={threads} -Xmx{params.mem_g}g \
            {params.extra} 2> {log}
        """
