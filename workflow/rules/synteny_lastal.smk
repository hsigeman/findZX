rule lastdb:
    input:
        synteny_ref
    output:
        outdir + "lastdb/" + "lastdb_" + synteny_abbr + ".prj"
    params:
        db_name= "lastdb_" + synteny_abbr,
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
        ref = ref_genome,
        db_ext = outdir + "lastdb/" + "lastdb_" + synteny_abbr + ".prj",
    output:
        outdir + "synteny_lastal/" + synteny_abbr + "_align"
    params:
        db = outdir + "lastdb/" + "lastdb_" + synteny_abbr
    threads: threads_max
    shell:
        """
        lastal -P {threads} {params.db} {input.ref}  | last-split > {output}
        """

rule maf_convert_syns:
    input:
        outdir + "synteny_lastal/" + synteny_abbr + "_align"
    output:
        outdir + "synteny_lastal/" + synteny_abbr + "_align_converted"
    threads: 1
    shell:
        """
        maf-convert psl {input} > {output}
        """