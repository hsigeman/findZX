rule lastal_syns:
    input:
        ref_genome
    output:
        outdir + "synteny_lastal/ + "_align"
    params:
        db = SYNS_DB
    threads: 15
    shell:
        """
        lastal {params.db} {input} -P 15 | last-split > {output}
        """

rule maf_convert_syns:
    input:
        COMP_GEN_SYNS + REF_NAME + "_align"
    output:
        COMP_GEN_SYNS + REF_NAME + "_align_converted"
    threads: 1
    shell:
        """
        maf-convert psl {input} > {output}
        """