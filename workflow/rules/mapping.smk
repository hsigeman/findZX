rule map_reads:
    input:
        reads= get_trimmed_reads if config['trim_reads'] | config['trim_and_subsample'] | config['subsample_only'] else get_fastq_new,
        idx=rules.bwa_index.output,
    output:
        temp(map_dir + "{sample}__{group}.sorted.tmp.bam"),
    log:
        logs_dir + "bwa_mem/{sample}__{group}.log",
    message:
        "Mapping reads to reference genome. Sample: {wildcards.sample}"
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
        map_dir + "{sample}__{group}.sorted.tmp.bam",
    output:
        temp(map_dir + "{sample}__{group}.sorted.bam"),
    log:
        logs_dir + "samtools_view/{sample}__{group}.log",
    message:
        "Filtering mapped reads, only proper pairs and quality score >20. Sample: {wildcards.sample}"
    params:
        extra="-bf 0x2 -F 260 -q 20" # optional params string
    wrapper:
        "0.78.0/bio/samtools/view"


rule mark_duplicates:
    input:
        map_dir + "{sample}__{group}.sorted.bam",
    output:
        bam=dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.unfiltered.bam",
        metrics=qc_dir + "dedup/{sample}__{group}.metrics.txt",
    log:
        logs_dir + "picard/dedup/{sample}__{group}.log",
    message:
        "Remove duplicate reads. Sample: {wildcards.sample}"
    params:
        extra="REMOVE_DUPLICATES=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 TMP_DIR=" + dedup_dir + "/"
    resources:
        mem_mb=mem_max
    wrapper:
        "0.74.0/bio/picard/markduplicates"


rule bamtools_filter:
    input:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.unfiltered.bam",
    output:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.0.{ED, [0-9]+}.bam", 
    params:
        tags = ["NM:<={ED}"]
    log:
        logs_dir + "bamtools/{sample}-{group}.{ED}.log",
    message:
        "Filter mapped reads for mismatches. Maximum {wildcards.ED} mismatches. Sample: {wildcards.sample}"
    wrapper:
        "0.74.0/bio/bamtools/filter"


rule samtools_index:
    input:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.bam", 
    output:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.bam.bai", 
    log:
        logs_dir + "samtools/{sample}-{group}.{ED}.log",
    message:
        "Index BAM file: {wildcards.sample}__{wildcards.group}.sorted.dedup.mismatch.{wildcards.ED}.bam"
    wrapper:
        "0.74.0/bio/samtools/index"


rule samtools_stats:
    input:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.bam",
    output:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.samtools.stats.txt",
    params:
        extra="",                       # Optional: extra arguments.
    log:
        logs_dir + "samtools_stats/{sample}__{group}.{ED}.log",
    message:
        "Calculating BAM file statistics: {wildcards.sample}__{wildcards.group}.sorted.dedup.mismatch.{wildcards.ED}.bam"
    wrapper:
        "0.74.0/bio/samtools/stats"


rule calc_cov:
    input:
        fai = ref_genome + ".fai",
        stats = dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.samtools.stats.txt",
    output:
        report(dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}_mean_coverage.txt", category = "Coverage")
    shell:
        "workflow/scripts/calc_cov.sh {input.stats} {input.fai} > {output}"


