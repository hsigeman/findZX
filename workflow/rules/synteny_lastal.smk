rule lastdb:
    input:
        synteny_ref
    output:
#        lastdb_dir + "lastdb_" + synteny_ref_name_simple + ".prj",
        log = lastdb_dir + "lastdb_" + synteny_ref_name_simple + ".log",
    params:
        db_name= "lastdb_" + synteny_ref_name_simple,
        synteny_dir = lastdb_dir
    threads: threads_max
    conda: 
        "../envs/last.yaml"
    shell:
        """
        lastdb -cR11 -P {threads} {params.db_name} {input}
        mv {params.db_name}* {params.synteny_dir}
#        touch {output}
        echo "DONE" > {output.log}
        """

rule lastal_syns:
    input:
        ref = ref_genome,
        log = lastdb_dir + "lastdb_" + synteny_ref_name_simple + ".log",
#        db_ext = lastdb_dir + "lastdb_" + synteny_ref_name_simple + ".prj",
    output:
        temp(windowCalc_het + synteny_abbr + "_align")
    params:
        db = lastdb_dir + "lastdb_" + synteny_ref_name_simple
    threads: threads_max
    conda: 
        "../envs/last.yaml"
    shell:
        """
        parallel-fasta "lastal -P {threads} {params.db}  | last-split" < {input.ref} > {output}
        """


rule maf_convert_syns:
    input:
        windowCalc_het + synteny_abbr + "_align"
    output:
        windowCalc_het + synteny_abbr + "_align_converted"
    conda: 
        "../envs/last.yaml"
    shell:
        """
        maf-convert psl {input} > {output}
        """