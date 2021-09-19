rule map_reads:
    input:
        reads= get_trimmed_reads if config['trim_reads'] else get_fastq_new,
        idx=rules.bwa_index.output,
    output:
        temp(map_dir + "{sample}__{unit}.sorted.tmp.bam"),
    log:
        logs_dir + "bwa_mem/{sample}__{unit}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra=get_read_group,
        sort="samtools",
        sort_order="coordinate",
        sort_extra="-@ " + str(threads_max),
    threads: threads_max
    wrapper:
        "0.74.0/bio/bwa/mem"


rule samtools_view:
    input:
        map_dir + "{sample}__{unit}.sorted.tmp.bam",
    output:
        temp(map_dir + "{sample}__{unit}.sorted.bam"),
    log:
        logs_dir + "samtools_view/{sample}__{unit}.log",
    params:
        extra="-bf 0x2" # optional params string
    wrapper:
        "0.78.0/bio/samtools/view"


rule mark_duplicates:
    input:
        map_dir + "{sample}__{unit}.sorted.bam",
    output:
        bam=dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.unfiltered.bam",
        metrics=qc_dir + "dedup/{sample}__{unit}.metrics.txt",
    log:
        logs_dir + "picard/dedup/{sample}__{unit}.log",
    params:
        extra="REMOVE_DUPLICATES=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true TMP_DIR=" + dedup_dir + "/"
    resources:
        mem_mb=mem_max
    wrapper:
        "0.74.0/bio/picard/markduplicates"


rule bamtools_filter:
    input:
        dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.unfiltered.bam",
    output:
        dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.0.{ED, [0-9]+}.bam", 
    params:
        tags = ["NM:<={ED}"]
    log:
        logs_dir + "bamtools/{sample}-{unit}.{ED}.log",
    wrapper:
        "0.74.0/bio/bamtools/filter"


rule samtools_index:
    input:
        dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.{ED}.bam", 
    output:
        dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.{ED}.bam.bai", 
    log:
        logs_dir + "samtools/{sample}-{unit}.{ED}.log",
    wrapper:
        "0.74.0/bio/samtools/index"


rule samtools_stats:
    input:
        dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.{ED}.bam",
    output:
        dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.{ED}.samtools.stats.txt",
    params:
        extra="",                       # Optional: extra arguments.
    log:
        logs_dir + "samtools_stats/{sample}__{unit}.{ED}.log",
    wrapper:
        "0.74.0/bio/samtools/stats"


rule calc_cov:
    input:
        fai = ref_genome + ".fai",
        stats = dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.{ED}.samtools.stats.txt",
    output:
        report(dedup_dir + "{sample}__{unit}.sorted.dedup.mismatch.{ED}_mean_coverage.txt", category = "Coverage")
    shell:
        "code/calc_cov.sh {input.stats} {input.fai} > {output}"


