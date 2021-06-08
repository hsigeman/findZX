rule confirm_sexing:
    input:
        gencov = outdir + "coverage/" + "gencov.nodup.nm.{ED}.out",
        het = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.NR.bed",
        stats = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.{{ED}}.samtools.stats.txt", zip, u=units.itertuples())
    output:
        read_length = outdir + "output/no_synteny/plots/" + ".misc/" + "read_length.sorted.nodup.nm.{ED}.csv",
        gencov_het = report(outdir + "output/no_synteny/plots/" + "confirm_sexing_indv.{ED}.pdf", category="2. Confirm sexing", caption="../report/confirm_sexing.rst")
    threads: 1
    params:
        map_dir = dedup_dir + "*.sorted.dedup.nm.{ED}.samtools.stats.txt",
        hetero = expand("{u.sample}__{u.unit}", u=heterogametic.itertuples()),
        homo = expand("{u.sample}__{u.unit}", u=homogametic.itertuples()),
        chromosomes = CHROMOSOMES
    threads: 1
    conda: 
        "../envs/R.yaml"
    shell:
        """
        python code/read_length.py <(for FILE in $(ls {params.map_dir}); do echo \"${{FILE##*/}}\"; grep \"average length\" $FILE; done) > {output.read_length}

        Rscript code/histogram_indv.R {input.gencov} {input.het} {output.read_length} {output.gencov_het} no-synteny {params.chromosomes} {params.hetero} {params.homo}
        """


rule plotting:
    input:
        cov = expand(outdir + "output/no_synteny/tables/" + "gencov.nodup.nm.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.{bp}bp.out"
    output:
        out_scatter = report(outdir + "output/no_synteny/plots/scatter.{bp}bp.pdf", category="3. Output plots", caption="../report/scatter_plots.rst",),
        out = touch(outdir + "output/no_synteny/plots/" + ".misc/" + "plotting.{bp}bp.done")
    params:
        chromosomes = CHROMOSOMES,
	    chromosomes_highlight = CHROMOSOMES_HIGHLIGHT,
	    ED = expand("{ED}", ED = EDIT_DIST),
    conda: 
        "../envs/R.yaml"
    shell:
        """
        Rscript code/plot_windows.R {input.cov} {input.snp} {output.out_scatter} {params.chromosomes} {params.chromosomes_highlight} {params.ED} 
        """


rule plotting_linear:
    input:
        cov = expand(outdir + "output/no_synteny/tables/" + "gencov.nodup.nm.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.{bp}bp.out"
    output:
        absolute_out = report(outdir + "output/no_synteny/plots/sexSpecificValues.{bp}bp.pdf", category="3. Output plots", caption="../report/linear_plots.rst",),
        diff_out = report(outdir + "output/no_synteny/plots/sexDifference.{bp}bp.pdf", category="3. Output plots", caption="../report/linear_plots.rst",),
        out = touch(outdir + "output/no_synteny/plots/" + ".misc/" + "plotting.linear.{bp}bp.done")
    params:
        chromosomes = CHROMOSOMES,
	    ED = expand("{ED}", ED = EDIT_DIST),
	    nr_chromosomes = 50
    conda: 
        "../envs/R.yaml"
    shell:
        """
        Rscript code/plot_windows_linear.R {input.cov} {input.snp} {output.absolute_out} {output.diff_out} {params.chromosomes} {params.ED} {params.nr_chromosomes}
        """


rule plotting_chr:
    input:
        cov = expand(outdir + "output/no_synteny/tables/" + "gencov.nodup.nm.{ED}.chr.out", ED = EDIT_DIST),
        snp = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.chr.out"
    output:
        out_scatter2D = report(outdir + "output/no_synteny/plots/chr_scatter2D.pdf", category="3. Output plots", caption="../report/scatter_plots.rst",),
        out_scatter3D = report(outdir + "output/no_synteny/plots/chr_scatter3D.pdf", category="3. Output plots", caption="../report/scatter_plots.rst",),
        out = touch(outdir + "output/no_synteny/plots/" + ".misc/" + "plotting_chr.done")
    params:
        chromosomes = CHROMOSOMES,
	    ED = expand("{ED}", ED = EDIT_DIST)
    conda: 
        "../envs/R.yaml"
    shell:
        """
        Rscript code/scatterplot_chr.R {input.cov} {input.snp} {output.out_scatter2D} {output.out_scatter3D} {params.chromosomes} {params.ED}
        """
