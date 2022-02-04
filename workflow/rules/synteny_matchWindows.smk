rule matchScaffold2Chr:
    input:
        syns = windowCalc_het + synteny_abbr + "_align_converted",
        gencov = cov_dir + "genome_5kb_windows.out"
    output:
        windows = windowCalc_het + "genome_windows.out",
        bestMatch = windowCalc_het + "bestMatch.list",
        log = windowCalc_het + "bestMatch.status"
    params:
        abswindow = dir_path + "/" + windowCalc_het + "genome_windows.out",
        absBestMatch = dir_path + "/" + windowCalc_het + "bestMatch.unfilter.list",
        okWindows = dir_path + "/" + windowCalc_het + "okWindows.list",
        absBestMatchFilter = dir_path + "/" + windowCalc_het + "bestMatch.list",
        windowsfile = "genome_windows.out",
        absLog = dir_path + "/" + windowCalc_het + "bestMatch.status"
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
        cat {input.syns} | awk '{{print $10,$12,$13,$14,$16,$17,$1}}' | sed 's/ /\t/g' | bedtools intersect -a stdin -b {input.gencov} -wa -wb | awk '{{if($10-$9==\"5000\") print $8,$9,$10,$7,$1,$2,$3,$4,$5,$6}}' | sed 's/ /\t/g' | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' > {output.windows}

        sort -r -n -k2 < {output.windows} | awk '!x[$1]++' | sort -k1 | sed 's/STARTCOORD/\t/' | sed 's/ENDCOORD/\t/' > {params.absBestMatch}

        cat {params.absBestMatch} | cut -f 8-10 | sort | uniq -c | awk '$1<=2 {{print}}' | awk '{{print $2,$3,$4}}' | sed 's/ /\t/g' | sort | uniq | bedtools sort > {params.okWindows}

        bedtools intersect -a <(<{params.absBestMatch} awk '{{print $8,$9,$10,$0}}' | sed 's/ /\t/g' | sort -k1,1 | bedtools sort) -b {params.okWindows} -f 1 -r -wa | cut -f 4- | awk '$4>500 {{print}}'> {params.absBestMatchFilter}

        echo "DONE" > {params.absLog}
        """


rule matchScaffold2Chr_snp:
    input:
        bestMatch = windowCalc_het + "bestMatch.list",
        het = vcf_dir + ref_genome_name_simple + ".heterozygosity.5kb.windows.NR.bed",
        het_sexAverage = vcf_dir + ref_genome_name_simple + ".heterozygosity.sexAverage.NR.bed"
    output:
        bestMatch = windowCalc_het + "heterozygosity.bestMatch",
        bestMatch_small = windowCalc_het + "heterozygosity.bestMatch.small",
        bestMatch_sexAverage = windowCalc_het + "heterozygosity.bestMatch.sexAverage.bed",
        bestMatch_small_sexAverage = windowCalc_het + "heterozygosity.bestMatch.small.sexAverage.bed"
    params:
        hetero = expand("het:{u.sample}__{u.group}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.group}", u=homogametic.itertuples())
    threads: 1
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.bestMatch} -b {input.het} -wa -wb > {output.bestMatch}

        cut -f 8,9,10,14- {output.bestMatch} > {output.bestMatch_small}
        
        bedtools intersect -a {input.bestMatch} -b {input.het_sexAverage} -wa -wb > {output.bestMatch_sexAverage}

        cut -f 8,9,10,14- {output.bestMatch_sexAverage} > {output.bestMatch_small_sexAverage}
        """


rule matchScaffold2Chr_cov:
    input:
        bestMatch = windowCalc_het + "bestMatch.list",
        cov = cov_dir + "gencov.mismatch.{ED}.out",
        cov_sexAverage = cov_dir + "gencov.mismatch.{ED}.norm.sexAverage.out"
    output:
        bestMatch = windowCalc_het + "gencov.mismatch.{ED}.out",
        bestMatch_small = windowCalc_het + "gencov.mismatch.{ED}.small.out",
        bestMatch_sexAverage = windowCalc_het + "gencov.mismatch.{ED}.norm.sexAverage.out",
        bestMatch_small_sexAverage = windowCalc_het + "gencov.mismatch.{ED}.norm.sexAverage.small.out"
    params:
        hetero = expand("het:{u.sample}__{u.group}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.group}", u=homogametic.itertuples())
    threads: 1
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.bestMatch} -b {input.cov} -wa -wb > {output.bestMatch}
        cut -f 8,9,10,14- {output.bestMatch} > {output.bestMatch_small}
        bedtools intersect -a {input.bestMatch} -b {input.cov_sexAverage} -wa -wb > {output.bestMatch_sexAverage}
        cut -f 8,9,10,14- {output.bestMatch_sexAverage} > {output.bestMatch_small_sexAverage}
        """

rule synteny_stats:
    input:
        bestMatch = windowCalc_het + "bestMatch.list",
        ref_stats = qc_dir + "assembly_stats/" + ref_genome_name_simple + "_stats.txt",
    output:
        report(windowCalc_het + "synteny_stats.out", category="01 Sample and reference genome statistics", caption="../report/synteny_perc.rst"),
    threads: 1
    shell:
        """
        ref_length=$(cat {input.ref_stats} | cut -f 2 | tail -n 1) 
        match_bp=$(cat {input.bestMatch} | cut -f 1-3 | sort | uniq | awk '{{print $3-$2}}' | paste -sd+ - | bc )
        echo "scale=2; ${{match_bp}}/${{ref_length}}" | bc | awk '{{printf "%f", $0}}'> {output}
        """