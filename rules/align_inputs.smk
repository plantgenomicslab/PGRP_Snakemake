# Input routing for align_*/alignRSEM_* rules.
#
# When config["bbduk_enable"] is true, alignment consumes BBDuk clean output;
# otherwise it consumes fastp trim output directly. Helpers here keep that
# branching out of Snakefile so the four align rules just declare
# `input: unpack(align_inputs)`.
#
# Expects LAYOUT and PAIR_LIST to be defined by the parent Snakefile before
# this file is included.

BBDUK_ENABLE = config.get("bbduk_enable", False)

if BBDUK_ENABLE:
    include: "bbduk_filter.smk"


def _trim_inputs(wc):
    if LAYOUT == "PAIRED":
        return {
            "fwd_fastq": f"output/{wc.replicate}/{wc.sample}/trim/{wc.sample}{PAIR_LIST[0]}_trimmed.fq.gz",
            "rev_fastq": f"output/{wc.replicate}/{wc.sample}/trim/{wc.sample}{PAIR_LIST[1]}_trimmed.fq.gz",
        }
    return {"fwd_fastq": f"output/{wc.replicate}/{wc.sample}/trim/{wc.sample}_trimmed.fq.gz"}


def _bbduk_inputs(wc):
    if LAYOUT == "PAIRED":
        return {
            "fwd_fastq": f"output/{wc.replicate}/{wc.sample}/bbduk/{wc.sample}{PAIR_LIST[0]}_clean.fq.gz",
            "rev_fastq": f"output/{wc.replicate}/{wc.sample}/bbduk/{wc.sample}{PAIR_LIST[1]}_clean.fq.gz",
        }
    return {"fwd_fastq": f"output/{wc.replicate}/{wc.sample}/bbduk/{wc.sample}_clean.fq.gz"}


def align_inputs(wc):
    return _bbduk_inputs(wc) if BBDUK_ENABLE else _trim_inputs(wc)
