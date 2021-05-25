rule genome_faidx:
    input:
        ref_genome,
    output:
        ref_genome + ".fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "0.74.0/bio/samtools/faidx"

rule bwa_index:
    input:
        ref_genome,
    output:
        multiext(ref_genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "0.74.0/bio/bwa/index"