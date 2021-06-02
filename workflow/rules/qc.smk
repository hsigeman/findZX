rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html=qc_dir + "fastqc/{sample}__{unit}.untrimmed.html",
        zip=qc_dir + "fastqc/{sample}__{unit}.untrimmed_fastqc.zip"
    log:
        logs_dir + "fastqc/{sample}__{unit}.untrimmed.log",
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc:
    input:
        expand(qc_dir + "fastqc/{u.sample}__{u.unit}.untrimmed_fastqc.zip", u=units.itertuples()),
    output:
        report(qc_dir + "fastqc/multiqc.untrimmed.html"),
    log:
        logs_dir + "fastqc/multiqc.log",
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
    wrapper:
        "0.74.0/bio/fastqc"


checkpoint multiqc_2:
    input:
        expand(qc_dir + "fastqc/{u.sample}__{u.unit}.trimmed_fastqc.zip",
            u=units.itertuples()),
    output:
        report(qc_dir + "fastqc/multiqc.trimmed.html"),
    log:
        logs_dir + "fastqc/multiqc.trimmed.log",
    wrapper:
        "0.74.0/bio/multiqc"
