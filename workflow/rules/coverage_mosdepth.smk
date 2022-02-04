rule gencov_prepare_fasta:
    input:
        ref_genome + ".fai",
    output:
        filter_fai = cov_dir + ref_genome_name_simple + ".filter." + MIN_SIZE_SCAFFOLD + ".fasta.fai",
        windows = cov_dir + "genome_5kb_windows.out"
    params:
        MIN_SIZE_SCAFFOLD
    conda: 
        "../envs/bedtools.yaml"
    message:
        "Make 5kb genome windows for BEDtools"
    shell:
        """
	    cat {input} | awk '$2>= {params} {{print $0}}' > {output.filter_fai}
        bedtools makewindows -g {output.filter_fai} -w 5000 -s 5000 | awk '$3-$2==5000 {{print}}' | bedtools sort > {output.windows}
        """


rule mosdepth_make_bed:
    input:
        ref_genome + ".fai",
    output:
        filter_bed = cov_dir + ref_genome_name_simple + ".filter." + MIN_SIZE_SCAFFOLD + ".fasta.bed"
    params:
        MIN_SIZE_SCAFFOLD
    message:
        "Make filtered BED file"
    shell:
        """
        cat {input} | awk -v OFS="\t" '$2>= {params} {{print $1,"0",$2}}' > {output.filter_bed}
        """


min_cov=str(config['min_cov'])
max_cov=str(config['max_cov'])


rule mosdepth_by_threshold:
    input:
        bam=dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.bam",
        bai=dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.bam.bai",
        windows=cov_dir + "genome_5kb_windows.out",
    output:
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.mosdepth.global.dist.txt",
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.mosdepth.region.dist.txt",
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.regions.bed.gz",
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.thresholds.bed.gz",
        summary=cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.mosdepth.summary.txt", 
    log:
        logs_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.log",
    conda: 
        "../envs/mosdepth.yaml"
    threads: 3
    params:
        thresholds= min_cov + "," + max_cov, 
        prefix=cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}"
    shell:
        """
        mosdepth -t {threads} -b {input.windows} -n -x -Q 20 -T {params.thresholds} {params.prefix} {input.bam}
        """


rule mosdepth_plot:
    input:
        expand(cov_dir + "mosdepth_by_threshold/{u.sample}__{u.group}.mismatch.{{ED}}.mosdepth.global.dist.txt", u=units.itertuples()),
    output:
        out_all=cov_dir + "mosdepth_by_threshold/perSampleCoverage_{ED}.html",
    conda: 
        "../envs/python_gawk.yaml"
    shell:
        """
        python workflow/scripts/plot-dist.py {input} -o {output} 
        """

rule bedops_merge_prep: 
    input:
        cov=cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.regions.bed.gz",
        thresholds=cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.thresholds.bed.gz",
    output:
        cov=temp(cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.regions.bed"),
        thresholds=temp(cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.thresholds.bed"),
        cov_filtered_tmp=cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.regions.filtered.tmp",
        cov_filtered=cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.regions.filtered.bed"
    conda: 
        "../envs/bedtools.yaml"
    params: 
        min_cov_perc=0.9,
        max_cov_perc=0.1,
    shell:
        """
        gunzip -c {input.cov} | sort -k1,1 -k2,2g > {output.cov}
        gunzip -c {input.thresholds} | sort -k1,1 -k2,2g | awk -v min={params.min_cov_perc} '$5/5000>min {{print}}' | awk -v max={params.max_cov_perc} '$6/5000<max {{print}}'> {output.thresholds}
        bedtools intersect -a {output.cov} -b {output.thresholds} -wa > {output.cov_filtered_tmp}
        bedtools intersect -a {output.cov} -b {output.thresholds} -v -wa | awk '{{print $1,$2,$3,"NaN"}}' | sed 's/ /\t/g' >> {output.cov_filtered_tmp}
        cat {output.cov_filtered_tmp} | bedtools sort | awk '$3-$2==5000 {{print}}' > {output.cov_filtered}
        """

rule merge_bedfiles: 
    input:
        bed_hetero=expand(cov_dir + "mosdepth_by_threshold/{u.sample}__{u.group}.mismatch.{{ED}}.regions.filtered.bed",zip, u=heterogametic.itertuples()),
        bed_homo=expand(cov_dir + "mosdepth_by_threshold/{u.sample}__{u.group}.mismatch.{{ED}}.regions.filtered.bed",zip, u=homogametic.itertuples()),
    output:
        out=cov_dir + "gencov.mismatch.{ED}.out"
    params:
        min_cov=min_cov,
        max_cov=max_cov
    conda: 
        "../envs/bedtools.yaml"
    shell:
        """
        bedtools unionbedg -i {input.bed_hetero} {input.bed_homo} | awk '{{for(i=4;i<=NF;i++)if($i>{params.max_cov})$i="NaN"}}1' | awk '{{for(i=4;i<=NF;i++)if($i<{params.min_cov})$i="NaN"}}1' | sed 's/ /\t/g' > {output.out}
        """

rule normalize_cov_mean:
    input:
        cov_dir + "gencov.mismatch.{ED}.out"
    output:
        cov_dir + "gencov.mismatch.{ED}.norm.sexAverage.out"
    params:
        hetero = expand("het:{u.sample}__{u.group}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.group}", u=homogametic.itertuples())
    conda: 
        "../envs/python_gawk.yaml"
    message:
        "Normalize genome coverage values between samples, and calculate mean value per sex for each 5 kb window"
    shell:
        """
        python3 workflow/scripts/normalize_genCov.py {input} no-synteny {params.hetero} {params.homo} | awk '{{print $1,int($2),int($3),$4,$5}}' | sed 's/ /\t/g' > {output}
        """
