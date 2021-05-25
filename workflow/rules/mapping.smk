
rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        temp(outdir + "mapped/{sample}-{unit}.sorted.bam"),
    log:
        logs_dir + "bwa_mem/{sample}-{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
    threads: 2
    wrapper:
        "0.74.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        outdir + "mapped/{sample}-{unit}.sorted.bam",
    output:
        bam=outdir + "dedup/{sample}-{unit}.bam",
        metrics=qc_dir + "dedup/{sample}-{unit}.metrics.txt",
    log:
        logs_dir + "picard/dedup/{sample}-{unit}.log",
    params:
        "REMOVE_DUPLICATES=true"
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/markduplicates"

rule samtools_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/samtools/index/{prefix}.log",
    wrapper:
        "0.74.0/bio/samtools/index"