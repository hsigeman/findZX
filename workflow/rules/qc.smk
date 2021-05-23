rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html="results/qc/fastqc/{sample}-{heterogamety}.html",
        zip="results/qc/fastqc/{sample}-{heterogamety}.zip",
    log:
        "logs/fastqc/{sample}-{heterogamety}.log",
    wrapper:
        "0.74.0/bio/fastqc"


rule samtools_stats:
    input:
        "results/recal/{sample}-{heterogamety}.bam",
    output:
        "results/qc/samtools-stats/{sample}-{heterogamety}.txt",
    log:
        "logs/samtools-stats/{sample}-{heterogamety}.log",
    wrapper:
        "0.74.0/bio/samtools/stats"


rule multiqc:
    input:
        expand(
            [
                "results/qc/samtools-stats/{u.sample}-{u.heterogamety}.txt",
                "results/qc/fastqc/{u.sample}-{u.heterogamety}.zip",
                "results/qc/dedup/{u.sample}-{u.heterogamety}.metrics.txt",
            ],
            u=units.itertuples(),
        ),
    output:
        report(
            "results/qc/multiqc.html",
            caption="../report/multiqc.rst",
            category="Quality control",
        ),
    log:
        "logs/multiqc.log",
    wrapper:
        "0.74.0/bio/multiqc"
