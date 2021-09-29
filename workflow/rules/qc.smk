rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html=qc_dir + "fastqc/{sample}__{unit}.untrimmed.html",
        zip=qc_dir + "fastqc/{sample}__{unit}.untrimmed_fastqc.zip"
    log:
        logs_dir + "fastqc/{sample}__{unit}.untrimmed.log",
    message:
        "Running FastQC on the untrimmed fastq file. Sample: {wildcards.sample}"
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc:
    input:
        expand(qc_dir + "fastqc/{u.sample}__{u.unit}.untrimmed_fastqc.zip", u=units.itertuples()),
    output:
        report(qc_dir + "fastqc/multiqc.untrimmed.html", category="1. Quality control", caption="../report/multiqc_untrimmed.rst",),
    log:
        logs_dir + "fastqc/multiqc.log",
    message:
        "Running MultiQC on the untrimmed fastq files"
    wrapper:
        "0.74.0/bio/multiqc"


rule fastqc_2:
    input:
        unpack(get_trimmed_reads),
    output:
        html=qc_dir + "fastqc/{sample}__{unit}.trimmed.html",
        zip=qc_dir + "fastqc/{sample}__{unit}.trimmed_fastqc.zip"
    log:
        logs_dir + "fastqc/{sample}__{unit}.trimmed.log",
    message:
        "Running FastQC on the trimmed fastq file. Sample: {wildcards.sample}"
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc_2:
    input:
        expand(qc_dir + "fastqc/{u.sample}__{u.unit}.trimmed_fastqc.zip",
            u=units.itertuples()),
    output:
        report(qc_dir + "fastqc/multiqc.trimmed.html", category="1. Quality control", caption="../report/multiqc_trimmed.rst",),
    log:
        logs_dir + "fastqc/multiqc.trimmed.log",
    message:
        "Running MultiQC on the untrimmed fastq files"
    wrapper:
        "0.74.0/bio/multiqc"
