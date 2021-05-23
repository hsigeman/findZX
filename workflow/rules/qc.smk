rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html="results/qc/fastqc/{sample}-{unit}.html",
        zip="results/qc/fastqc/{sample}-{unit}.zip",
    log:
        "logs/fastqc/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}-{unit}.zip")
    output:
        "results/qc/multiqc.html"
    log:
        "logs/multiqc.log",
    wrapper:
        "0.74.0/bio/multiqc"
