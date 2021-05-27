rule confirm_sexing:
    input:
        gencov = outdir + "synteny_lastal/" + "gencov.nodup.nm.{ED}.small.out",
        het = outdir + "synteny_lastal/" + "heterozygosity.bestMatch.small",
        stats = expand(qc_dir + "dedup/{u.sample}__{u.unit}.sorted.dedup.nm.{{ED}}.samtools.stats.txt", zip, u=units.itertuples())
    output:
        read_length = outdir + "output/synteny/" + synteny_abbr + "/.misc/" + "read_length.sorted.nodup.nm.{ED}.csv",
        gencov_het = outdir + "output/synteny/" + synteny_abbr + "/plots/" + "confirm_sexing_indv.{ED}.pdf"
    threads: 1
    params:
        map_dir = qc_dir + "dedup/" + "*.sorted.dedup.nm.{ED}.samtools.stats.txt",
        hetero = expand("{u.sample}__{u.unit}", u=heterogametic.itertuples()),
        homo = expand("{u.sample}__{u.unit}", u=homogametic.itertuples()),
        chromosomes = CHROMOSOMES
    threads: 1
    shell:
        """
        python code/read_length.py <(for FILE in $(ls {params.map_dir}); do echo \"${{FILE##*/}}\"; grep \"average length\" $FILE; done) > {output.read_length}

        Rscript code/histogram_indv.R {input.gencov} {input.het} {output.read_length} {output.gencov_het} synteny {params.chromosomes} {params.hetero} {params.homo}
        """


rule plotting:
    input:
        cov = expand(outdir + "output/synteny/" + synteny_abbr + "/tables/" + "gencov.nodup.nm.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "diffHeterozygosity.{bp}bp.out"
    output:
        touch(outdir + "output/synteny/" + synteny_abbr + "/.misc/" +  "plotting.{bp}bp.done")
    params:
        out_scatter = outdir + "output/synteny/" + synteny_abbr + "/plots/scatter.{bp}bp.pdf",
        chromosomes = CHROMOSOMES,
	    chromosomes_highlight = CHROMOSOMES_HIGHLIGHT,
	    ED = expand("{ED}", ED = EDIT_DIST),
    shell:
        """
        Rscript code/plot_windows.R {input.cov} {input.snp} {params.out_scatter} {params.chromosomes} {params.chromosomes_highlight} {params.ED} 
        """


rule plotting_linear:
    input:
        cov = expand(outdir + "output/synteny/" + synteny_abbr + "/tables/" + "gencov.nodup.nm.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "diffHeterozygosity.{bp}bp.out"
    output:
        touch(outdir + "output/synteny/" + synteny_abbr + "/.misc/" + "plotting.linear.{bp}bp.done")
    params:
        absolute_out = outdir + "output/synteny/" + synteny_abbr + "/plots/sexSpecificValues.{bp}bp.pdf",
        diff_out = outdir + "output/synteny/" + synteny_abbr + "/plots/sexDifference.{bp}bp.pdf",
        chromosomes = CHROMOSOMES,
	    ED = expand("{ED}", ED = EDIT_DIST),
	    nr_chromosomes = 50
    shell:
        """
        Rscript code/plot_windows_linear.R {input.cov} {input.snp} {params.absolute_out} {params.diff_out} {params.chromosomes} {params.ED} {params.nr_chromosomes}
        """


rule plotting_chr:
    input:
        cov = expand(outdir + "output/synteny/" + synteny_abbr + "/tables/" + "gencov.nodup.nm.{ED}.chr.out", ED = EDIT_DIST),
        snp = outdir + "output/synteny/" + synteny_abbr + "/tables/" + "diffHeterozygosity.chr.out"
    output:
        touch(outdir + "output/synteny/" + synteny_abbr + "/.misc/" + "plotting_chr.done")
    params:
        out_scatter2D = protected(outdir + "output/synteny/" + synteny_abbr + "/plots/chr_scatter2D.pdf"),
        out_scatter3D = protected(outdir + "output/synteny/" + synteny_abbr + "/plots/chr_scatter3D.pdf"),
        chromosomes = CHROMOSOMES,
	    ED = expand("{ED}", ED = EDIT_DIST)
    shell:
        """
        Rscript code/scatterplot_chr.R {input.cov} {input.snp} {params.out_scatter2D} {params.out_scatter3D} {params.chromosomes} {params.ED}
        """
