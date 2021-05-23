
rule trimmomatic_pe:
    input:
        r1=RAW_FQ_DIR + "{forward}",
        r2=RAW_FQ_DIR + "{reverse}"
    output:
        r1=TRIMMED_FQ_DIR + "{S}.1.fastq.gz",
        r2=TRIMMED_FQ_DIR + "{S}.2.fastq.gz",
        # reads where trimming entirely removed the mate
        r1_unpaired=TRIMMED_FQ_DIR + "{S}.1.unpaired.fastq.gz",
        r2_unpaired=TRIMMED_FQ_DIR + "{S}.2.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["TRAILING:3"],
        # optional parameters
        extra="",
        compression_level="-9"
    threads:
        10
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/trimmomatic/pe"
        
rule fastqc_pre_trimming:
    input:
        r1=RAW_FQ_DIR + "{S}.1.fastq.gz",
        r2=RAW_FQ_DIR + "{S}.2.fastq.gz"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 1
    wrapper:
        "0.74.0/bio/fastqc"        

        
rule fastqc_post_trimming:
    input:
        r1=TRIMMED_FQ_DIR + "{S}.1.fastq.gz",
        r2=TRIMMED_FQ_DIR + "{S}.2.fastq.gz"
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc/{sample}.log"
    threads: 1
    wrapper:
        "0.74.0/bio/fastqc"          
