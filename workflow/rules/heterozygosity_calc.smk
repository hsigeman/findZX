rule proportion_heterozygosity:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    output:
        het = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.bed")
    log: 
        outdir + "variant_calling/" + ref_genome_name_simple + ".proportion_heterozygosity.log"
    params:
        hetero = expand("het:{u.sample}__{u.unit}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.unit}", u=homogametic.itertuples())
    resources:
        mem_mb=2048
    shell:
        """
        python3 code/heterozygosity_per_indv.py {input} {output.het} {params.hetero} {params.homo} > {log}
        """


rule proportion_heterozygosity_window:
    input:
        het = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.bed",
        windows = outdir + "coverage/" + "genome_5kb_windows.out"
    output:
        het_sorted = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.sorted.bed"),
        windows_sorted = temp(outdir + "variant_calling/" + "genome_5kb_windows.out"),
        het_sorted_window = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.bed")
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools sort -i {input.windows} > {output.windows_sorted}
        bedtools sort -i {input.het} > {output.het_sorted}
        bedtools intersect -a {input.windows} -b {output.het_sorted} -wa -wb -sorted | cut -f 1-3,7- > {output.het_sorted_window}
        """


rule proportion_heterozygosity_window_2:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.bed"
    output:
        het_sorted_window_mean = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.NR.bed",
	    het_sexAverage = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.sexAverage.NR.bed"
    params:
        hetero = expand("het:{u.sample}__{u.unit}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.unit}", u=homogametic.itertuples())
    conda: 
        "../envs/python_gawk.yaml"
    shell:
        """
        code/sum_heterozygosity_per_5kb.sh {input} > {output.het_sorted_window_mean}
	    
        python3 code/mean_heterozygosity_per_sex.py {output.het_sorted_window_mean} no-synteny {params.hetero} {params.homo} > {output.het_sexAverage}
	    """
