rule map_reads:
    input:
        reads=get_trimmed_reads,
        idx=rules.bwa_index.output,
    output:
        temp(outdir + "mapped/{sample}__{unit}.sorted.bam"),
    log:
        logs_dir + "bwa_mem/{sample}__{unit}.log",
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
        outdir + "mapped/{sample}__{unit}.sorted.bam",
    output:
        bam=outdir + "dedup/{sample}__{unit}.sorted.dedup.nm.all.bam",
        metrics=qc_dir + "dedup/{sample}__{unit}.metrics.txt",
    log:
        logs_dir + "picard/dedup/{sample}__{unit}.log",
    params:
        extra="REMOVE_DUPLICATES=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/picard/markduplicates"


rule bamtools_filter:
    input:
        outdir + "dedup/{sample}__{unit}.sorted.dedup.nm.all.bam",
    output:
        outdir + "dedup/{sample}__{unit}.sorted.dedup.nm.0.{ED, [0-9]+}.bam", 
    params:
        tags = ["NM:<={ED}"]
    log:
        logs_dir + "bamtools/{sample}-{unit}.{ED}.log",
    wrapper:
        "0.74.0/bio/bamtools/filter"


rule samtools_index:
    input:
        outdir + "dedup/{sample}__{unit}.sorted.dedup.nm.{ED}.bam", 
    output:
        outdir + "dedup/{sample}__{unit}.sorted.dedup.nm.{ED}.bam.bai", 
    log:
        logs_dir + "samtools/{sample}-{unit}.{ED}.log",
    wrapper:
        "0.74.0/bio/samtools/index"


rule samtools_stats:
    input:
        outdir + "dedup/{sample}__{unit}.sorted.dedup.nm.{ED}.bam",
    output:
        qc_dir + "dedup/{sample}__{unit}.sorted.dedup.nm.{ED}.samtools.stats.txt",
    params:
        extra="",                       # Optional: extra arguments.
    log:
        logs_dir + "samtools_stats/{sample}__{unit}.{ED}.log",
    wrapper:
        "0.74.0/bio/samtools/stats"