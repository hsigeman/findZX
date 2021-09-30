checkpoint genome_faidx:
    input:
        ref_genome,
    output:
        ref_genome + ".fai",
    log:
        logs_dir + "samtools/genome-faidx.log",
    cache: True
    message:
        "Index reference genome (samtools)"
    wrapper:
        "0.74.0/bio/samtools/faidx"



rule MapCaller_index:
    input:
        ref_genome,
    output:
        ref_genome_dir + "/" + ref_genome_name_simple + ".bwt",
    log:
        logs_dir + "mapcaller/mapcaller_index.log",
    params: 
        name=ref_genome_name_simple,
        dir=ref_genome_dir + "/"
    conda: 
        "../envs/mapcaller.yaml"
    message:
        "Index reference genome (MapCaller)"
    shell:
        """
        MapCaller index {input} {params.name}
        mv {params.name}* {params.dir}
        touch {output}
        """
