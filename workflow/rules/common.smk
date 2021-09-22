import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version

min_version("5.18.0")



###### Config file and sample sheets #####


#configfile: "config/config.yml"


#samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

units = pd.read_table(config["units"], dtype=str).set_index(
    ["sample", "unit"], drop=False
)


units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels]
)  # enforce str in index


homogametic = units[units["unit"] == "homogametic"]
heterogametic = units[units["unit"] == "heterogametic"]


##### Wildcard constraints #####


wildcard_constraints:
    sample="|".join(units["sample"]),
    unit="|".join(units["unit"]),




##### Helper functions #####


# This function does not work with newer pandas versions. Update. 
def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {"r1": fastqs.fq1, "r2": fastqs.fq2}
    return {"r1": fastqs.fq1}


def get_contigs():
    with checkpoints.genome_faidx.get().output[0].open() as fai:
        return pd.read_table(fai, header=None, usecols=[0], squeeze=True, dtype=str)

def is_single_end(sample, unit):
    """Return True if sample-unit is single end."""
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_trimmed_reads(wildcards):
    """Get trimmed reads of given sample-unit."""
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand(
            outdir + "trimmed/{sample}__{unit}.{group}.fastq.gz",
            group=[1, 2],
            **wildcards
        )
    # single end sample
    return outdir + "trimmed/{sample}__{unit}.fastq.gz".format(**wildcards)


def get_fastq_new(wildcards):
    """Get fastq files of given sample-unit."""
    fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    if len(fastqs) == 2:
        return {fastqs.fq1, fastqs.fq2}
    return {fastqs.fq1}


def get_read_group(wildcards):
    """Denote sample name and platform in read group."""
    return r"-M -R '@RG\tID:{sample}__{unit}\tSM:{sample}__{unit}'".format(**wildcards)


def get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        outdir + "dedup/{sample}__{unit}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
    )


def new_get_sample_bams(wildcards):
    """Get all aligned reads of given sample."""
    return expand(
        outdir + "dedup/{sample}__{unit}.sorted.dedup.nm.{ED}.bam",
        sample=wildcards.sample,
        unit=units.loc[wildcards.sample].unit,
        ED = EDIT_DIST
    )

def get_fastq_r1(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1"]].dropna().values.flatten()

def get_fastq_r2(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq2"]].dropna().values.flatten()