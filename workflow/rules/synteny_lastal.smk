rule lastdb:
    input:
        synteny_ref
    output:
        outdir + "lastdb/" + "lastdb_" + synteny_ref_name_simple + ".prj"
    params:
        db_name= "lastdb_" + synteny_ref_name_simple,
        synteny_dir = outdir + "lastdb/"
    threads: 15
    shell:
        """
        lastdb -cR11 -P 15 {params.db_name} {input}
        mv {params.db_name}.* {params.synteny_dir}
        touch {output}
        """


rule lastal_syns:
    input:
        ref_genome
    output:
        outdir + "synteny_lastal/ + synteny_ref_name_simple + "_align"
    params:
        db = outdir + "lastdb/" + "lastdb_" + synteny_ref_name_simple
    threads: 15
    shell:
        """
        lastal {params.db} {input} -P 15 | last-split > {output}
        """

rule maf_convert_syns:
    input:
        outdir + "synteny_lastal/ + synteny_ref_name_simple + "_align"
    output:
        outdir + "synteny_lastal/ + synteny_ref_name_simple + "_align_converted"
    threads: 1
    shell:
        """
        maf-convert psl {input} > {output}
        """