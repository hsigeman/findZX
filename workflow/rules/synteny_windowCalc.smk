rule calculate_heterozygosity:
    input:
        outdir + "synteny_lastal/" + "heterozygosity.bestMatch.small.sexAverage.bed"
    output:
        Mb = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "diffHeterozygosity.1Mbp.out",
        kb = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "diffHeterozygosity.100kbp.out",
        chr = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "diffHeterozygosity.chr.out"
    threads: 1
#    params:
#        chromosomes = CHROMOSOMES
    shell:
        """
        Rscript code/calculate_windows.R {input} {output.Mb} {output.kb} {output.chr}
        """


rule calculate_ratio:
    input:
        outdir + "synteny_lastal/" + "gencov.nodup.nm.{ED}.norm.sexAverage.small.out",
    output:
        Mb = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "gencov.nodup.nm.{ED}.1Mbp.out",
        kb = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "gencov.nodup.nm.{ED}.100kbp.out",
        chr = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "gencov.nodup.nm.{ED}.chr.out"
    threads: 1
#    params:
#        chromosomes = CHROMOSOMES
    shell:
        """
        Rscript code/calculate_windows.R {input} {output.Mb} {output.kb} {output.chr}
        """
