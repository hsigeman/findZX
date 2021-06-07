rule gencov_prepare_fasta:
    input:
        ref_genome + ".fai",
    output:
        filter_fai = outdir + "coverage/" + ref_genome_name_simple + ".filter." + MIN_SIZE_SCAFFOLD + ".fasta.fai",
        windows = outdir + "coverage/" + "genome_5kb_windows.out"
    params:
        MIN_SIZE_SCAFFOLD
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
	    cat {input} | awk '$2>= {params} {{print $0}}' > {output.filter_fai}
        bedtools makewindows -g {output.filter_fai} -w 5000 -s 5000 | awk '$3-$2==5000 {{print}}' | bedtools sort > {output.windows}
        """


rule gencov_bedtools:
    input:
        bam_hetero = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.{{ED}}.bam", zip, u=heterogametic.itertuples()),
        bam_homo = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.{{ED}}.bam", zip, u=homogametic.itertuples()),
        bai_hetero = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.{{ED}}.bam.bai", zip, u=heterogametic.itertuples()),
        bai_homo = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.{{ED}}.bam.bai", zip, u=homogametic.itertuples()),
        bed = outdir + "coverage/" + "genome_5kb_windows.out"
    output:
        outdir + "coverage/" + "gencov.nodup.nm.{ED}.out"
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools multicov -bams {input.bam_hetero} {input.bam_homo} -bed {input.bed} -p -q 20 > {output}
        """


rule normalize_cov_mean:
    input:
        outdir + "coverage/" + "gencov.nodup.nm.{ED}.out"
    output:
        outdir + "coverage/" + "gencov.nodup.nm.{ED}.norm.sexAverage.out"
    params:
        hetero = expand("het:{u.sample}__{u.unit}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.unit}", u=homogametic.itertuples())
    conda: 
        "../envs/python_gawk.yaml"
    shell:
        """
        python3 code/normalize_genCov.py {input} no-synteny {params.hetero} {params.homo} > {output}
        """