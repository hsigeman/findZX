rule map_reads:
    input:
        log=logs_dir + "mapcaller/mapcaller_index.log",
        fq1=outdir + "trimmed/{sample}__{group}.1.fastq.gz" if config['trim_reads'] | config['trim_and_subsample'] | config['subsample_only'] else get_fastq_r1,
        fq2=outdir + "trimmed/{sample}__{group}.2.fastq.gz" if config['trim_reads'] | config['trim_and_subsample'] | config['subsample_only'] else get_fastq_r2,
    output:
        bam=temp(dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.unfiltered.unsorted.bam"),
        vcf=temp(vcf_dir + "{sample}__{group}.vcf")
    conda: 
        "../envs/mapcaller.yaml"
    params:
        ref=ref_genome_dir + "/" + ref_genome_name_simple + "_mapcaller",
        id="{sample}__{group}"
    log:
        logs_dir + "mapcaller/{sample}__{group}.log",
    threads: threads_max
    shell:
        "MapCaller -i {params.ref} -f {input.fq1} -f2 {input.fq2} -bam {output.bam} -vcf {output.vcf} -ad 1 -t {threads}"


rule add_readgroup:
    input:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.unfiltered.unsorted.bam"
    output:
        temp(dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.unfiltered.unsorted.RG.bam")
    params:
        get_read_group_MapCaller
    conda: 
        "../envs/samtools.yaml"
    log:
        logs_dir + "samtools_addRG/{sample}__{group}.log",
    message:
        "Add RG. Sample: {wildcards.sample}"
    shell:
        "samtools addreplacerg {input} -r {params} | samtools view -bS - > {output}"


rule samtools_sort:
    input:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.unfiltered.unsorted.RG.bam"
    output:
        dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.unfiltered.bam"
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads: threads_max
    log:
        logs_dir + "samtools_sort/{sample}__{group}.log",
    message:
        "Sorting reads. Sample: {wildcards.sample}"
    wrapper:
        "0.78.0/bio/samtools/sort"


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
    conda: 
        "../envs/samtools.yaml"
    log:
        logs_dir + "samtools/{sample}-{group}.{ED}.log",
    message:
        "Index BAM file: {wildcards.sample}__{wildcards.group}.sorted.dedup.mismatch.{wildcards.ED}.bam"
    shell:
        "samtools index {input} {output}"


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


