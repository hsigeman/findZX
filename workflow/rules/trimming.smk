if config['trim_reads']: 
    rule trim_reads_pe:
        input:
            unpack(get_fastq),
        output:
            r1=temp(outdir + "trimmed/{sample}__{unit}.1.fastq.gz"),
            r2=temp(outdir + "trimmed/{sample}__{unit}.2.fastq.gz"),
            r1_unpaired=temp(outdir + "trimmed/{sample}__{unit}.1.unpaired.fastq.gz"),
            r2_unpaired=temp(outdir + "trimmed/{sample}__{unit}.2.unpaired.fastq.gz"),
            trimlog=temp(outdir + "trimmed/{sample}__{unit}.trimlog.txt"),
        params:
            **config["params"]["trimmomatic"]["pe"],
            extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        threads:
            10
#        log:
#            logs_dir + "trimmomatic/{sample}__{unit}.log",
        wrapper:
            "0.74.0/bio/trimmomatic/pe"


if config['trim_reads']:
    rule trim_reads_se:
        input:
            unpack(get_fastq),
        output:
            temp(outdir + "trimmed/{sample}__{unit}.fastq.gz"),
        params:
            **config["params"]["trimmomatic"]["se"],
            extra="",
        threads:
            10
#        log:
#            logs_dir + "trimmomatic/{sample}__{unit}.log",
        wrapper:
            "0.74.0/bio/trimmomatic/se"


