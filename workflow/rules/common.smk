import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")



###### Config file and sample sheets #####


#configfile: "config/config.yml"


#samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "group"], drop=False
)


units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index


homogametic = units[units["group"] == "homogametic"]
heterogametic = units[units["group"] == "heterogametic"]


##### Wildcard constraints #####


wildcard_constraints:
    sample="|".join(units["sample"]),
    group="|".join(units["group"]),




##### Helper functions #####


# This function does not work with newer pandas versions. Update. 
def get_fastq(wildcards):
    """Get fastq files of given sample-group."""
    fastqs = units.loc[(wildcards.sample, wildcards.group), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)

def is_single_end(sample, group):
    """Return True if sample-group is single end."""
    return pd.isnull(units.loc[(sample, group), "fq2"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-group."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            outdir + "trimmed/{sample}__{group}.{nr}.fastq.gz",
            nr=[1, 2],
            **wildcards
        )
    # single end sample
    return outdir + "trimmed/{sample}__{group}.fastq.gz".format(**wildcards)


def get_fastq_new(wildcards):
    """Get fastq files of given sample-group."""
    fastqs = units.loc[(wildcards.sample, wildcards.group), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {fastqs.fq1, fastqs.fq2}
    return {fastqs.fq1}


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-M -R '@RG\tID:{sample}__{group}\tSM:{sample}__{group}'".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        outdir + "dedup/{sample}__{group}.bam",
        sample=wildcards.sample,
        group=units.loc[wildcards.sample].group,
    )


def new_get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        outdir + "dedup/{sample}__{group}.sorted.dedup.nm.{ED}.bam",
        sample=wildcards.sample,
        group=units.loc[wildcards.sample].group,
        ED = EDIT_DIST
    )

def get_fastq_r1(wildcards):
    return units.loc[(wildcards.sample, wildcards.group), ["fq1"]].dropna().values.flatten()

def get_fastq_r2(wildcards):
    return units.loc[(wildcards.sample, wildcards.group), ["fq2"]].dropna().values.flatten()