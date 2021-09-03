rule calculate_heterozygosity_chr:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.sexAverage.NR.bed"
    output:
        outdir + "output/no_synteny/tables/" + "diffHeterozygosity.chr.out"
    threads: 1
    conda: 
        "../envs/R.yaml"
    shell:
        """
        Rscript code/calculate_chr.R {input} {output}
        """


rule calculate_heterozygosity_window:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.sexAverage.NR.bed"
    output:
        outdir + "output/no_synteny/tables/" + "diffHeterozygosity.{window}bp.out"
    threads: 1
    params:
        "{window}"
    conda: 
        "../envs/R.yaml"
    shell:
        """
        Rscript code/calculate_windows_userSpec.R {input} {output} {params}
        """


rule calculate_ratio_chr:
    input:
        outdir + "coverage/" + "gencov.mismatch.{ED}.norm.sexAverage.out"
    output:
        outdir + "output/no_synteny/tables/" + "gencov.mismatch.{ED}.chr.out"
    threads: 1
    conda: 
        "../envs/R.yaml"
    shell:
        """
        Rscript code/calculate_chr.R {input} {output}
        """


rule calculate_ratio_window:
    input:
        outdir + "coverage/" + "gencov.mismatch.{ED}.norm.sexAverage.out"
    output:
        outdir + "output/no_synteny/tables/" + "gencov.mismatch.{ED}.{window}bp.out"
    threads: 1
    params:
        window="{window}"
    conda: 
        "../envs/R.yaml"
    shell:
        """
        Rscript code/calculate_windows_userSpec.R {input} {output} {params}
        """