rule proportion_heterozygosity:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    output:
        het = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.bed")
    log: 
        outdir + "variant_calling/" + ref_genome_name_simple + ".proportion_heterozygosity.log"
    params:
        hetero = expand("het:{u.sample}__{u.group}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.group}", u=homogametic.itertuples())
    message:
        "Scores individuals as heterozygotes or homozygotes for all sites in the VCF file"
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
        windows_sorted = temp(outdir + "variant_calling/" + "genome_5kb_windows.out"),
        het_tmp_window = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.tmp"),
        het_sorted_window = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.bed")
    params: 
        split = outdir + "variant_calling/split_file_",
        dir = outdir + "variant_calling"
    conda: 
        "../envs/bedtools.yaml"
    message:
        "Calculates number of heterozygous sites per individual and 5kb window"
    shell:
        """
        bedtools sort -i {input.windows} > {output.windows_sorted}
        split -l 100000 {input.het} {params.split}
        ls {params.dir} | grep split_file | while read file ; do bedtools intersect -a {input.windows} -b {params.dir}/${{file}} -wa -wb | cut -f 1-3,7- ; done > {output.het_tmp_window} 
        sort -k1,1 -k2,2 {output.het_tmp_window} > {output.het_sorted_window}
        rm {params.split}*
        """


rule proportion_heterozygosity_window_2:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.bed"
    output:
        het_sorted_window_mean = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.NR.bed",
        het_sexAverage = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.sexAverage.NR.bed"
    params:
        hetero = expand("het:{u.sample}__{u.group}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.group}", u=homogametic.itertuples())
    conda: 
        "../envs/python_gawk.yaml"
    message:
        "Calculates mean heterozygosity per sex and 5kb window"
    shell:
        """
        workflow/scripts/sum_heterozygosity_per_5kb.sh {input} > {output.het_sorted_window_mean}
	    
        python3 workflow/scripts/mean_heterozygosity_per_sex.py {output.het_sorted_window_mean} no-synteny {params.hetero} {params.homo} > {output.het_sexAverage}
        """
