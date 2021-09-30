if config['trim_reads']: 
    rule trim_reads_pe:
        input:
            unpack(get_fastq),
        output:
            r1=temp(outdir + "trimmed/{sample}__{group}.1.fastq.gz"),
            r2=temp(outdir + "trimmed/{sample}__{group}.2.fastq.gz"),
            r1_unpaired=temp(outdir + "trimmed/{sample}__{group}.1.unpaired.fastq.gz"),
            r2_unpaired=temp(outdir + "trimmed/{sample}__{group}.2.unpaired.fastq.gz"),
 #           trimlog=temp(outdir + "trimmed/{sample}__{group}.trimlog.txt"),
        params:
            **config["params"]["trimmomatic"]["pe"],
          #  extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        threads:
            10
        log:
            logs_dir + "trimmomatic/{sample}__{group}.log",
        message:
            "Trimming fastq reads with Trimmomatic. Sample: {wildcards.sample}"
        wrapper:
            "0.74.0/bio/trimmomatic/pe"


if config['trim_and_subsample']: 
    rule trim_reads_pe_2:
        input:
            unpack(get_fastq),
        output:
            r1=temp(outdir + "trimmed/{sample}__{group}.1.tmp.fastq.gz"),
            r2=temp(outdir + "trimmed/{sample}__{group}.2.tmp.fastq.gz"),
            r1_unpaired=temp(outdir + "trimmed/{sample}__{group}.1.unpaired.fastq.gz"),
            r2_unpaired=temp(outdir + "trimmed/{sample}__{group}.2.unpaired.fastq.gz"),
 #           trimlog=temp(outdir + "trimmed/{sample}__{group}.trimlog.txt"),
        params:
            **config["params"]["trimmomatic"]["pe"],
          #  extra=lambda w, output: "-trimlog {}".format(output.trimlog),
        threads:
            10
        log:
            logs_dir + "trimmomatic/{sample}__{group}.log",
        message:
            "Trimming fastq reads with Trimmomatic. Sample: {wildcards.sample}"
        wrapper:
            "0.74.0/bio/trimmomatic/pe"


if config['trim_and_subsample']: 
    rule subsample:
        input:
            fq1=outdir + "trimmed/{sample}__{group}.1.tmp.fastq.gz",
            fq2=outdir + "trimmed/{sample}__{group}.2.tmp.fastq.gz"
        output:
            r1=temp(outdir + "trimmed/{sample}__{group}.1.fastq.gz"),
            r2=temp(outdir + "trimmed/{sample}__{group}.2.fastq.gz"),
        params:
            basepairs= config["subsample_basepairs"],
        conda: 
            "../envs/bbtools.yaml"
        log: 
            logs_dir + "subsample/{sample}__{group}.log",
        message:
            "Subsample reads to {params.basepairs} bp. Sample: {wildcards.sample}"
        shell:
            """
            reformat.sh int=f in1={input.fq1} in2={input.fq2} out1={output.r1} out2={output.r2} samplebasestarget={params.basepairs} 2> {log}
            """ 

if config['subsample_only']: 
    rule subsample_2:
        input:
            fq1=get_fastq_r1,
            fq2=get_fastq_r2
        output:
            r1=temp(outdir + "trimmed/{sample}__{group}.1.fastq.gz"),
            r2=temp(outdir + "trimmed/{sample}__{group}.2.fastq.gz"),
        params:
            basepairs= config["subsample_basepairs"],
        conda: 
            "../envs/bbtools.yaml"
        log: 
            logs_dir + "subsample/{sample}__{group}.log",
        message:
            "Subsample reads to {params.basepairs} bp. Sample: {wildcards.sample}"
        shell:
            """
            reformat.sh int=f in1={input.fq1} in2={input.fq2} out1={output.r1} out2={output.r2} samplebasestarget={params.basepairs} 2> {log}
            """ 