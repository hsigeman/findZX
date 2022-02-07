num=range(1,999)

rule confirm_sexing:
    input:
        gencov = windowCalc_cov + "gencov.mismatch.{ED}.small.out",
        het = windowCalc_het + "heterozygosity.bestMatch.small",
        stats = expand(dedup_dir + "{u.sample}__{u.group}.sorted.dedup.mismatch.{{ED}}.samtools.stats.txt", zip, u=units.itertuples())
    output:
        read_length = plots_dir + ".misc/" + "read_length.sorted.nodup.mismatch.{ED}.csv",
        gencov_het = report(plots_dir + "6_confirmSexing.samplesSeparately.mismatch.{ED}.pdf", category="Output plot type 6 (\"confirm sexing\")", subcategory="mismatch option: {ED}", caption="../report/confirm_sexing.rst")
    threads: 1
    params:
        map_dir = dedup_dir + "*.sorted.dedup.mismatch.{ED}.samtools.stats.txt",
        hetero = expand("{u.sample}__{u.group}", u=heterogametic.itertuples()),
        homo = expand("{u.sample}__{u.group}", u=homogametic.itertuples()),
        chromosomes = CHROMOSOMES
    threads: 1
    conda: 
        "../envs/R.yaml"
    log:
        plot_log + "confirm_sexing_synteny.mismatch.{ED}.log"
    message:
        "Plotting results 6_confirmSexing"
    shell:
        """
        python workflow/scripts/read_length.py <(for FILE in $(ls {params.map_dir}); do echo \"${{FILE##*/}}\"; grep \"average length\" $FILE; done) > {output.read_length}

        Rscript workflow/scripts/histogram_indv.R {input.gencov} {input.het} {output.read_length} {output.gencov_het} synteny {params.chromosomes} {params.hetero} {params.homo} 2> {log}
        """


rule highlight_file:
    params: 
        highlight_chr=config['synteny_chr_highlight'] 
    output: 
        outputdir + "highlight_file.list"
    shell: 
        """
        echo {params.highlight_chr} | tr " " "\n" > {output}
        """


if not config['synteny_chr_highlight']: 
    rule plotting:
        input:
            cov = expand(tables_dir + "diffGenomeCoverage.mismatch.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
            snp = tables_dir + "diffHeterozygosity.{bp}bp.out",
            chromosomes_highlight = outputdir + "highlight_file.list"
        output:
            out_scatter = report(plots_dir + "4_sexDifferences.{bp}bp.pdf", category="Output plot type 4 (scatter plots)", subcategory = "window size: {bp} bp", caption="../report/scatter_plots.rst",),
            out_scatter_highlight = plots_dir + "4_sexDifferences.{bp}bp.highlight.pdf",
            out = touch(plots_dir + ".misc/" +  "plotting.{bp}bp.done")
        params:
            chromosomes = CHROMOSOMES,
            ED = expand("{ED}", ED = EDIT_DIST),
            window = "{bp}"
        conda: 
            "../envs/R.yaml"
        log:
            plot_log + "4_sexDifferences.mismatch.{bp}.log"
        message:
            "Plotting results 4_sexDifferences"
        shell:
            """
            Rscript workflow/scripts/plot_windows.R {input.cov} {input.snp} {output.out_scatter} {params.chromosomes} {input.chromosomes_highlight} {params.ED} {params.window} 2> {log}
            """

else:
    rule plotting:
        input:
            cov = expand(tables_dir + "diffGenomeCoverage.mismatch.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
            snp = tables_dir + "diffHeterozygosity.{bp}bp.out",
            chromosomes_highlight = outputdir + "highlight_file.list"
        output:
            out_scatter = report(plots_dir + "4_sexDifferences.{bp}bp.window.pdf", category="Output plot type 4 (scatter plots)", subcategory = "window size: {bp} bp", caption="../report/scatter_plots.rst",),
            out_scatter_highlight = report(plots_dir + "4_sexDifferences.{bp}bp.window.highlight.pdf", category="Output plot type 4 (scatter plots)", subcategory = "window size: {bp} bp", caption="../report/scatter_plots_highlight.rst"),
            out = touch(plots_dir + ".misc/" +  "plotting.{bp}bp.done")
        params:
            chromosomes = CHROMOSOMES,
            ED = expand("{ED}", ED = EDIT_DIST),
            window = "{bp}"
        conda: 
            "../envs/R.yaml"
        log:
            plot_log + "4_sexDifferences.mismatch.{bp}.log"
        message:
            "Plotting results 4_sexDifferences"
        shell:
            """
            Rscript workflow/scripts/plot_windows.R {input.cov} {input.snp} {output.out_scatter} {params.chromosomes} {input.chromosomes_highlight} {params.ED} {params.window} 2> {log}
            """


rule linear_models:
    input:
        cov = expand(tables_dir + "diffGenomeCoverage.mismatch.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = tables_dir + "diffHeterozygosity.{bp}bp.out",
        chromosomes_highlight = outputdir + "highlight_file.list"
    output:
        plot = report(plots_dir + "5_linear_model.plot.{bp}bp.pdf", category="Output plot type 5 (linear models)", subcategory = "window size: {bp} bp", caption="../report/linear_models.rst"),
        table = report(tables_dir + "linear_model_results_estimate_CI.{bp}bp.html", category = "Output tables", subcategory = "window size: {bp} bp", caption="../report/linear_models_html.rst"),
        out = touch(plots_dir + ".misc/" +  "linear.model.{bp}bp.done")
    params:
        chromosomes = CHROMOSOMES,
        nr_chromosomes = 50,
        ED = expand("{ED}", ED = EDIT_DIST),
        window = "{bp}"
    conda: 
        "../envs/R2.yaml"
    log:
        plot_log + "linear.models.{bp}.log"
    message:
        "Linear model analyses"
    shell:
        """
        Rscript workflow/scripts/linear_model_plotting_and_stats_test.R {input.cov} {input.snp} {output.plot} {output.table} {params.chromosomes} {input.chromosomes_highlight} {params.ED} {params.window} {params.nr_chromosomes} 2> {log}
        """

rule plotting_linear:
    input:
        cov = expand(tables_dir + "diffGenomeCoverage.mismatch.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = tables_dir + "diffHeterozygosity.{bp}bp.out"
    output:
        absolute_out = report(plots_dir + "2_sexesSeparate.genomeWide.{bp}bp.window.pdf", category="Output plot type 2 (genome-wide sexes separately)", subcategory = "window size: {bp} bp", caption="../report/linear_plots_sexesSeparate.rst",),
        diff_out = report(plots_dir + "1_sexDifferences.genomeWide.{bp}bp.window.pdf", category="Output plot type 1 (genome-wide sex differences)", subcategory = "window size: {bp} bp", caption="../report/linear_plots.rst",),
        out = touch(plots_dir + ".misc/" + "plotting.linear.{bp}bp.done")
    params:
        chromosomes = CHROMOSOMES,
	    ED = expand("{ED}", ED = EDIT_DIST),
	    nr_chromosomes = 50,
        window = "{bp}"
    conda: 
        "../envs/R.yaml"
    log:
        plot_log + "plotting_linear.{bp}bp.log"
    message:
        "Plotting results 1_sexDifferences"
    shell:
        """
        Rscript workflow/scripts/plot_windows_linear.R {input.cov} {input.snp} {output.absolute_out} {output.diff_out} {params.chromosomes} {params.ED} {params.nr_chromosomes} {params.window} 2> {log}
        """


rule plotting_chr:
    input:
        cov = expand(tables_dir + "diffGenomeCoverage.mismatch.{ED}.chr.out", ED = EDIT_DIST),
        snp = tables_dir + "diffHeterozygosity.chr.out"
    output:
        out_scatter2D = report(plots_dir + "3_sexDifferences.chromosome.pdf", category="Output plot type 3 (scatter plots with chromosome/scaffold length)", caption="../report/scatter_2D.rst",),
        out = touch(plots_dir + ".misc/" + "plotting_chr.done")
    params:
        chromosomes = CHROMOSOMES,
	    ED = expand("{ED}", ED = EDIT_DIST)
    conda: 
        "../envs/R.yaml"
    log:
        plot_log + "plotting_chr.log"
    message:
        "Plotting results 3_sexDifferences"
    shell:
        """
        Rscript workflow/scripts/scatterplot_chr.R {input.cov} {input.snp} {output.out_scatter2D} {params.chromosomes} {params.ED} 2> {log}
        """


rule table_readme:
    input:
        "workflow/report/output_table_README.md"
    output:
        tables_dir + "output_table_README.md"
    shell:
        """
        cp {input} {output}
        """