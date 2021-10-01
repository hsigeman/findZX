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
    message:
        "Make 5kb genome windows for BEDtools"
    shell:
        """
	    cat {input} | awk '$2>= {params} {{print $0}}' > {output.filter_fai}
        bedtools makewindows -g {output.filter_fai} -w 5000 -s 5000 | awk '$3-$2==5000 {{print}}' | bedtools sort > {output.windows}
        """


rule gencov_bedtools:
    input:
        bam_hetero = expand(dedup_dir + "{u.sample}__{u.group}.mapcaller.sorted.dedup.mismatch.{{ED}}.bam", zip, u=heterogametic.itertuples()),
        bam_homo = expand(dedup_dir + "{u.sample}__{u.group}.mapcaller.sorted.dedup.mismatch.{{ED}}.bam", zip, u=homogametic.itertuples()),
        bai_hetero = expand(dedup_dir + "{u.sample}__{u.group}.mapcaller.sorted.dedup.mismatch.{{ED}}.bam.bai", zip, u=heterogametic.itertuples()),
        bai_homo = expand(dedup_dir + "{u.sample}__{u.group}.mapcaller.sorted.dedup.mismatch.{{ED}}.bam.bai", zip, u=homogametic.itertuples()),
        bed = outdir + "coverage/" + "genome_5kb_windows.out"
    output:
        outdir + "coverage/" + "gencov.mismatch.{ED}.out"
    conda: 
        "../envs/bedtools.yaml"
    message:
        "Calculate genome coverage with BEDtools"
    shell:
        """
        bedtools multicov -bams {input.bam_hetero} {input.bam_homo} -bed {input.bed} > {output}
        """


rule normalize_cov_mean:
    input:
        outdir + "coverage/" + "gencov.mismatch.{ED}.out"
    output:
        outdir + "coverage/" + "gencov.mismatch.{ED}.norm.sexAverage.out"
    params:
        hetero = expand("het:{u.sample}__{u.group}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.group}", u=homogametic.itertuples())
    conda: 
        "../envs/python_gawk.yaml"
    message:
        "Normalize genome coverage values between samples, and calculate mean value per sex for each 5 kb window"
    shell:
        """
        python3 workflow/scripts/normalize_genCov.py {input} no-synteny {params.hetero} {params.homo} > {output}
        """