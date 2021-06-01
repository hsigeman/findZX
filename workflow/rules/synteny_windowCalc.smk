rule calculate_heterozygosity_chr:
    input:
        outdir + "synteny_lastal/" + "heterozygosity.bestMatch.small.sexAverage.bed"
    output:
        outdir + "output/synteny/" + synteny_abbr + "/tables/" + "diffHeterozygosity.chr.out"
    threads: 1
    shell:
        """
        Rscript code/calculate_chr.R {input} {output}
        """


rule calculate_heterozygosity_window:
    input:
        outdir + "synteny_lastal/" + "heterozygosity.bestMatch.small.sexAverage.bed"
    output:
        outdir + "output/synteny/" + synteny_abbr + "/tables/" + "diffHeterozygosity.{window}bp.out"
    threads: 1
    params:
        "{window}"
    shell:
        """
        Rscript code/calculate_windows_userSpec.R {input} {output} {params}
        """


rule calculate_ratio_chr:
    input:
        outdir + "synteny_lastal/" + "gencov.nodup.nm.{ED}.norm.sexAverage.small.out",
    output:
        outdir + "output/synteny/" + synteny_abbr + "/tables/" + "gencov.nodup.nm.{ED}.chr.out"
    threads: 1
    shell:
        """
        Rscript code/calculate_chr.R {input} {output}
        """


rule calculate_ratio_window:
    input:
        outdir + "synteny_lastal/" + "gencov.nodup.nm.{ED}.norm.sexAverage.small.out",
    output:
        outdir + "output/synteny/" + synteny_abbr + "/tables/" + "gencov.nodup.nm.{ED}.{window}bp.out"
    threads: 1
    params:
        window="{window}"
    shell:
        """
        Rscript code/calculate_windows_userSpec.R {input} {output} {params}
        """