rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}.zip",
    log:
        "logs/fastqc/{sample}.log",
    wrapper:
        "0.74.0/bio/fastqc"

rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}.zip", sample=samples['sample'])
    output:
 #       report(
            "results/qc/multiqc.html"
 #           caption="../report/multiqc.rst",
  #          category="Quality control",
  #      ),
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "0.74.0/bio/multiqc"