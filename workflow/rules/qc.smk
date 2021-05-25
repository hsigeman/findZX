rule fastqc:
    input:
        unpack(get_fastq),
    output:
        html=qc_dir + "fastqc/{sample}-{unit}.untrimmed.html",
        zip=qc_dir + "fastqc/{sample}-{unit}.untrimmed_fastqc.zip"
    log:
        logs_dir + "fastqc/{sample}-{unit}.untrimmed.log",
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc:
    input:
        expand(qc_dir + "fastqc/{u.sample}-{u.unit}.untrimmed_fastqc.zip", u=units.itertuples()),
    output:
        report(qc_dir + "fastqc/multiqc.untrimmed.html"),
    log:
        logs_dir + "fastqc/multiqc.log",
    wrapper:
        "0.74.0/bio/multiqc"


rule trim_reads_pe:
    input:
        unpack(get_fastq),
    output:
        r1=temp(outdir + "trimmed/{sample}-{unit}.1.fastq.gz"),
        r2=temp(outdir + "trimmed/{sample}-{unit}.2.fastq.gz"),
        r1_unpaired=temp(outdir + "trimmed/{sample}-{unit}.1.unpaired.fastq.gz"),
        r2_unpaired=temp(outdir + "trimmed/{sample}-{unit}.2.unpaired.fastq.gz"),
        trimlog=outdir + "trimmed/{sample}-{unit}.trimlog.txt",
    params:
        **config["params"]["trimmomatic"]["pe"],
        extra=lambda w, output: "-trimlog {}".format(output.trimlog),
    log:
        logs_dir + "trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/trimmomatic/pe"


rule trim_reads_se:
    input:
        unpack(get_fastq),
    output:
        temp(outdir + "trimmed/{sample}-{unit}.fastq.gz"),
    params:
        **config["params"]["trimmomatic"]["se"],
        extra="",
    log:
        logs_dir + "trimmomatic/{sample}-{unit}.log",
    wrapper:
        "0.74.0/bio/trimmomatic/se"


rule fastqc_2:
    input:
        unpack(get_trimmed_reads),
    output:
        html=qc_dir + "fastqc/{sample}-{unit}.trimmed.html",
        zip=qc_dir + "fastqc/{sample}-{unit}.trimmed_fastqc.zip"
    log:
        logs_dir + "fastqc/{sample}-{unit}.trimmed.log",
    wrapper:
        "0.74.0/bio/fastqc"


rule multiqc_2:
    input:
        expand(qc_dir + "fastqc/{u.sample}-{u.unit}.trimmed_fastqc.zip",
            u=units.itertuples()),
    output:
        report(qc_dir + "fastqc/multiqc.trimmed.html"),
    log:
        logs_dir + "fastqc/multiqc.trimmed.log",
    wrapper:
        "0.74.0/bio/multiqc"
