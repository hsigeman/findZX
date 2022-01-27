rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html=qc_dir + "fastqc/{sample}__{group}.untrimmed.html",
        zip=qc_dir + "fastqc/{sample}__{group}.untrimmed_fastqc.zip"
    log:
        logs_dir + "fastqc/{sample}__{group}.untrimmed.log",
    message:
        "Running FastQC on the untrimmed fastq file. Sample: {wildcards.sample}"
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc:
    input:
        expand(qc_dir + "fastqc/{u.sample}__{u.group}.untrimmed_fastqc.zip", u=units.itertuples()),
    output:
        report(qc_dir + "fastqc/multiqc.untrimmed.html", category="1. Quality control", caption="../report/multiqc_untrimmed.rst",),
    log:
        logs_dir + "fastqc/multiqc.log",
    message:
        "Running MultiQC on the untrimmed fastq files"
    params: 
        qc_dir + "fastqc/multiqc.untrimmed.html"
    conda: 
        "../envs/multiqc.yaml"
    shell:
        "multiqc {input} -n {params} -f"


rule fastqc_2:
    input:
        unpack(get_trimmed_reads),
    output:
        html=qc_dir + "fastqc/{sample}__{group}.trimmed.html",
        zip=qc_dir + "fastqc/{sample}__{group}.trimmed_fastqc.zip"
    log:
        logs_dir + "fastqc/{sample}__{group}.trimmed.log",
    message:
        "Running FastQC on the trimmed fastq file. Sample: {wildcards.sample}"
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc_2:
    input:
        expand(qc_dir + "fastqc/{u.sample}__{u.group}.trimmed_fastqc.zip",
            u=units.itertuples()),
    output:
        report(qc_dir + "fastqc/multiqc.trimmed.html", category="1. Quality control", caption="../report/multiqc_trimmed.rst",),
    log:
        logs_dir + "fastqc/multiqc.trimmed.log",
    params: 
        qc_dir + "fastqc/multiqc.trimmed.html"
    message:
        "Running MultiQC on the trimmed fastq files"
    conda: 
        "../envs/multiqc.yaml"
    shell:
        "multiqc {input} -n {params} -f"


rule run_assembly_stats:
    input:
        assembly=ref_genome,
    output:
        assembly_stats=qc_dir + "assembly_stats/" + ref_genome_name_simple + "_stats.txt",
    params:
        extra="-s -l " + "MIN_SIZE_SCAFFOLD",
    log:
        logs_dir + "assembly_stats/" + ref_genome_name_simple + ".assembly-stats.log",
    threads: 1
    wrapper:
        "v0.86.0/bio/assembly-stats"