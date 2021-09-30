rule matchScaffold2Chr:
    input:
        syns = outdir + "synteny_lastal/" + synteny_abbr + "/" + synteny_abbr + "_align_converted",
        gencov = outdir + "coverage/" + "genome_5kb_windows.out"
    output:
        windows = outdir + "synteny_lastal/" + synteny_abbr + "/"  + "genome_windows.out",
        bestMatch = outdir + "synteny_lastal/" + synteny_abbr + "/" + "bestMatch.list",
        log = outdir + "synteny_lastal/" + synteny_abbr + "/" + "bestMatch.status"
    params:
        abswindow = dir_path + "/" + outdir + "synteny_lastal/" + synteny_abbr + "/" + "genome_windows.out",
        absBestMatch = dir_path + "/" + outdir + "synteny_lastal/" + synteny_abbr + "/" + "bestMatch.unfilter.list",
        okWindows = dir_path + "/" + outdir + "synteny_lastal/" + synteny_abbr + "/" + "okWindows.list",
        absBestMatchFilter = dir_path + "/" + outdir + "synteny_lastal/" + synteny_abbr + "/" + "bestMatch.list",
        windowsfile = "genome_windows.out",
        absLog = dir_path + "/" + outdir + "synteny_lastal/" + synteny_abbr + "/" + "bestMatch.status"
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
        bestMatch = outdir + "synteny_lastal/" + synteny_abbr + "/" + "bestMatch.list",
        het = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.NR.bed",
        het_sexAverage = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.sexAverage.NR.bed"
    output:
        bestMatch = outdir + "synteny_lastal/" + synteny_abbr + "/" + "heterozygosity.bestMatch",
        bestMatch_small = outdir + "synteny_lastal/" + synteny_abbr + "/" + "heterozygosity.bestMatch.small",
        bestMatch_sexAverage = outdir + "synteny_lastal/" + synteny_abbr + "/" + "heterozygosity.bestMatch.sexAverage.bed",
        bestMatch_small_sexAverage = outdir + "synteny_lastal/" + synteny_abbr + "/" + "heterozygosity.bestMatch.small.sexAverage.bed"
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
        bestMatch = outdir + "synteny_lastal/" + synteny_abbr + "/" + "bestMatch.list",
        cov = outdir + "coverage/" + "gencov.mismatch.{ED}.out",
        cov_sexAverage = outdir + "coverage/" + "gencov.mismatch.{ED}.norm.sexAverage.out"
    output:
        bestMatch = outdir + "synteny_lastal/" + synteny_abbr + "/" + "gencov.mismatch.{ED}.out",
        bestMatch_small = outdir + "synteny_lastal/" + synteny_abbr + "/" + "gencov.mismatch.{ED}.small.out",
        bestMatch_sexAverage = outdir + "synteny_lastal/" + synteny_abbr + "/" + "gencov.mismatch.{ED}.norm.sexAverage.out",
        bestMatch_small_sexAverage = outdir + "synteny_lastal/" + synteny_abbr + "/" + "gencov.mismatch.{ED}.norm.sexAverage.small.out"
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