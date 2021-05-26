##########################################################
################ STATISTICAL CALCULATIONS ################
##########################################################

rule calculate_heterozygosity:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.sexAverage.NR.bed"
    output:
        Mb = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.1Mbp.out",
        kb = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.100kbp.out",
        chr = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.chr.out"
    threads: 1
#    params:
#        chromosomes = CHROMOSOMES
    shell:
        """
        Rscript code/calculate_windows.R {input} {output.Mb} {output.kb} {output.chr}
        """

rule calculate_ratio:
    input:
        outdir + "coverage/" + "gencov.nodup.nm.{ED}.norm.sexAverage.out"
    output:
        Mb = outdir + "output/no_synteny/tables/" + "gencov.nodup.nm.{ED}.1Mbp.out",
        kb = outdir + "output/no_synteny/tables/" + "gencov.nodup.nm.{ED}.100kbp.out",
        chr = outdir + "output/no_synteny/tables/" + "gencov.nodup.nm.{ED}.chr.out"
    threads: 1
#    params:
#        chromosomes = CHROMOSOMES
    shell:
        """
        Rscript code/calculate_windows.R {input} {output.Mb} {output.kb} {output.chr}
        """
