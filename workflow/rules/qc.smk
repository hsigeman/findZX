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

#rule fastqc_unzip:
#    input:
#        "results/qc/fastqc/{sample}-{unit}.zip"
#    output:
#       log="results/qc/fastqc/{sample}-{unit}.done.txt",
#        directory("results/qc/fastqc/{sample}-{unit}")
#    shell:
#        """
#        unzip -o -d {output} {input}
#        touch {output.log}
#        """


rule multiqc:
    input:
        expand(
            [
                "results/qc/fastqc/{u.sample}-{u.unit}.html",
            ],
            u=units.itertuples(),
        ),
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