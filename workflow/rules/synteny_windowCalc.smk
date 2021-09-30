rule calculate_heterozygosity_chr:
    input:
        windowCalc_het + "heterozygosity.bestMatch.small.sexAverage.bed"
    output:
        tables_dir + "diffHeterozygosity.chr.out"
    threads: 1
    conda: 
        "../envs/R.yaml"
    message: 
        "Calculating heterozygosity sex differences per chromosome."
    shell:
        """
        Rscript workflow/scripts/calculate_chr.R {input} {output}
        """


rule calculate_heterozygosity_window:
    input:
        windowCalc_het + "heterozygosity.bestMatch.small.sexAverage.bed"
    output:
        tables_dir + "diffHeterozygosity.{window}bp.out"
    threads: 1
    params:
        "{window}"
    conda: 
        "../envs/R.yaml"
    message: 
        "Calculating heterozygosity sex differences. Window size: {wildcards.window}"
    shell:
        """
        Rscript workflow/scripts/calculate_windows_userSpec.R {input} {output} {params}
        """


rule calculate_ratio_chr:
    input:
        windowCalc_het + "gencov.mismatch.{ED}.norm.sexAverage.small.out",
    output:
        tables_dir + "diffGenomeCoverage.mismatch.{ED}.chr.out"
    threads: 1
    conda: 
        "../envs/R.yaml"
    message: 
        "Calculating genome coverage sex differences. Mismatch setting: {wildcards.ED}. Per chromosome."
    shell:
        """
        Rscript workflow/scripts/calculate_chr.R {input} {output}
        """


rule calculate_ratio_window:
    input:
        windowCalc_cov + "gencov.mismatch.{ED}.norm.sexAverage.small.out",
    output:
        tables_dir + "diffGenomeCoverage.mismatch.{ED}.{window}bp.out"
    threads: 1
    conda: 
        "../envs/R.yaml"
    params:
        window="{window}"
    message: 
        "Calculating genome coverage sex differences. Mismatch setting: {wildcards.ED}. Window size: {wildcards.window}"
    shell:
        """
        Rscript workflow/scripts/calculate_windows_userSpec.R {input} {output} {params}
        """