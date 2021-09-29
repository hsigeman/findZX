rule confirm_sexing:
    input:
        gencov = outdir + "coverage/" + "gencov.mismatch.{ED}.out",
        het = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.NR.bed",
        stats = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.mismatch.{{ED}}.samtools.stats.txt", zip, u=units.itertuples())
    output:
        read_length = outdir + "output/no_synteny/plots/" + ".misc/" + "read_length.sorted.nodup.mismatch.{ED}.csv",
        gencov_het = outdir + "output/no_synteny/plots/" + "5_confirmSexing.samplesSeparately.mismatch.{ED}.pdf"
    threads: 1
    params:
        map_dir = dedup_dir + "*.sorted.dedup.mismatch.{ED}.samtools.stats.txt",
        hetero = expand("{u.sample}__{u.unit}", u=heterogametic.itertuples()),
        homo = expand("{u.sample}__{u.unit}", u=homogametic.itertuples()),
        chromosomes = CHROMOSOMES
    threads: 1
    conda: 
        "../envs/R.yaml"
    log:
        logs_dir + "plotting/confirm_sexing.mismatch.{ED}.log"
    message:
        "Plotting results 5_confirmSexing"
    shell:
        """
        python code/read_length.py <(for FILE in $(ls {params.map_dir}); do echo \"${{FILE##*/}}\"; grep \"average length\" $FILE; done) > {output.read_length}

        Rscript code/histogram_indv.R {input.gencov} {input.het} {output.read_length} {output.gencov_het} no-synteny {params.chromosomes} {params.hetero} {params.homo} 2> {log}
        """


rule highlight_file:
    params: 
        highlight_chr=config['chr_highlight'] 
    output: 
        outdir + "output/no_synteny/" + "highlight_file.list"
    shell: 
        """
        echo {params.highlight_chr} | tr " " "\n" > {output}
        """


if not config['chr_highlight']:
    rule plotting:
        input:
            cov = expand(outdir + "output/no_synteny/tables/" + "diffGenomeCoverage.mismatch.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
            snp = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.{bp}bp.out",
            chromosomes_highlight = outdir + "output/no_synteny/" + "highlight_file.list",
        output:
            out_scatter = report(outdir + "output/no_synteny/plots/4_sexDifferences.{bp}bp.pdf", category="3. Output plots", caption="../report/scatter_plots.rst",),
            out_scatter_highlight = outdir + "output/no_synteny/plots/4_sexDifferences.{bp}bp.highlight.pdf",
            out = touch(outdir + "output/no_synteny/plots/" + ".misc/" + "plotting.{bp}bp.done")
        params:
            chromosomes = CHROMOSOMES,
	        ED = expand("{ED}", ED = EDIT_DIST),
            window = "{bp}"
        conda: 
            "../envs/R.yaml"
        log:
            logs_dir + "plotting/4_sexDifferences.mismatch.{bp}.log"
        message:
            "Plotting results 4_sexDifferences"
        shell:
            """
            Rscript code/plot_windows.R {input.cov} {input.snp} {output.out_scatter} {params.chromosomes} {input.chromosomes_highlight} {params.ED} {params.window} 2> {log}
            """

else: 
    rule plotting:
        input:
            cov = expand(outdir + "output/no_synteny/tables/" + "diffGenomeCoverage.mismatch.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
            snp = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.{bp}bp.out",
            chromosomes_highlight = outdir + "output/no_synteny/" + "highlight_file.list",
        output:
            out_scatter = report(outdir + "output/no_synteny/plots/4_sexDifferences.{bp}bp.window.pdf", category="3. Output plots", caption="../report/scatter_plots.rst",),
            out_scatter_highlight = report(outdir + "output/no_synteny/plots/4_sexDifferences.{bp}bp.window.highlight.pdf", category="3. Output plots", caption="../report/scatter_plots_highlight.rst"),
            out = touch(outdir + "output/no_synteny/plots/" + ".misc/" + "plotting.{bp}bp.done")
        params:
            chromosomes = CHROMOSOMES,
            ED = expand("{ED}", ED = EDIT_DIST),
            window = "{bp}"
        conda: 
            "../envs/R.yaml"
        log:
            logs_dir + "plotting/4_sexDifferences.mismatch.{bp}.log"
        message:
            "Plotting results 4_sexDifferences"
        shell:
            """
            Rscript code/plot_windows.R {input.cov} {input.snp} {output.out_scatter} {params.chromosomes} {input.chromosomes_highlight} {params.ED} {params.window} 2> {log}
            """


rule plotting_linear:
    input:
        cov = expand(outdir + "output/no_synteny/tables/" + "diffGenomeCoverage.mismatch.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.{bp}bp.out"
    output:
        absolute_out = report(outdir + "output/no_synteny/plots/2_sexesSeparate.genomeWide.{bp}bp.window.pdf", category="3. Output plots", caption="../report/linear_plots.rst",),
        diff_out = report(outdir + "output/no_synteny/plots/1_sexDifferences.genomeWide.{bp}bp.window.pdf", category="3. Output plots", caption="../report/linear_plots.rst",),
        out = touch(outdir + "output/no_synteny/plots/" + ".misc/" + "plotting.linear.{bp}bp.done")
    params:
        chromosomes = CHROMOSOMES,
	    ED = expand("{ED}", ED = EDIT_DIST),
	    nr_chromosomes = 50,
        window = "{bp}"
    conda: 
        "../envs/R.yaml"
    log:
        logs_dir + "plotting/plotting_linear.{bp}bp.log"
    message:
        "Plotting results 1_sexDifferences"
    shell:
        """
        Rscript code/plot_windows_linear.R {input.cov} {input.snp} {output.absolute_out} {output.diff_out} {params.chromosomes} {params.ED} {params.nr_chromosomes} {params.window} 2> {log}
        """


rule plotting_chr:
    input:
        cov = expand(outdir + "output/no_synteny/tables/" + "diffGenomeCoverage.mismatch.{ED}.chr.out", ED = EDIT_DIST),
        snp = outdir + "output/no_synteny/tables/" + "diffHeterozygosity.chr.out"
    output:
        out_scatter2D = report(outdir + "output/no_synteny/plots/3_sexDifferences.chromosome.pdf", category="3. Output plots", caption="../report/scatter_2D.rst",),
 #       out_scatter3D = report(outdir + "output/no_synteny/plots/chr_scatter3D.pdf", category="3. Output plots", caption="../report/scatter_3D.rst",),
        out = touch(outdir + "output/no_synteny/plots/" + ".misc/" + "plotting_chr.done")
    params:
        chromosomes = CHROMOSOMES,
	    ED = expand("{ED}", ED = EDIT_DIST)
    conda: 
        "../envs/R.yaml"
    log:
        logs_dir + "plotting/plotting_chr.log"
    message:
        "Plotting results 3_sexDifferences"
    shell:
        """
        Rscript code/scatterplot_chr.R {input.cov} {input.snp} {output.out_scatter2D} {params.chromosomes} {params.ED} 2> {log}
        """


rule table_readme:
    input:
        "config/output_table_README.md"
    output:
        outdir + "output/no_synteny/tables/" + "output_table_README.md"
    shell:
        """
        cp {input} {output}
        """