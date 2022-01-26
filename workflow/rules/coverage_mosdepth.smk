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


rule mosdepth_by_threshold:
    input:
        bam=dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.bam",
        bai=dedup_dir + "{sample}__{group}.sorted.dedup.mismatch.{ED}.bam.bai",
    output:
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.mosdepth.global.dist.txt",
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.mosdepth.region.dist.txt",
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.regions.bed.gz",
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.thresholds.bed.gz",
        summary=cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.mosdepth.summary.txt", 
    log:
        logs_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.log",
    params:
        by="5000",  # optional, window size,  specifies --by for mosdepth.region.dist.txt and regions.bed.gz
        thresholds="1,5,10,30",  # optional, specifies --thresholds for thresholds.bed.gz
        extra= "--no-per-base -Q 20 -x --use-median"
    # additional decompression threads through `--threads`
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v0.86.0/bio/mosdepth"


rule mosdepth_plot:
    input:
        expand(cov_dir + "mosdepth_by_threshold/{u.sample}__{u.group}.mismatch.{{ED}}.mosdepth.global.dist.txt", u=units.itertuples()),
    output:
        cov_dir + "mosdepth_by_threshold/perSampleCoverage_{ED}.html",
    conda: 
        "../envs/python_gawk.yaml"
    shell:
        "python workflow/scripts/plot-dist.py {input} -o {output}" 


rule bedops_merge_prep: 
    input:
        cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.regions.bed.gz"
    output:
        temp(cov_dir + "mosdepth_by_threshold/{sample}__{group}.mismatch.{ED}.regions.bed")
    shell:
        "gunzip -c {input} | sort -k1,1 -k2,2g > {output}"


rule bedops_merge: 
    input:
        bed_hetero=expand(cov_dir + "mosdepth_by_threshold/{u.sample}__{u.group}.mismatch.{{ED}}.regions.bed",zip, u=heterogametic.itertuples()),
        bed_homo=expand(cov_dir + "mosdepth_by_threshold/{u.sample}__{u.group}.mismatch.{{ED}}.regions.bed",zip, u=homogametic.itertuples()),
        windows=cov_dir + "genome_5kb_windows.out"
    output:
        union=temp(cov_dir + "mosdepth_by_threshold/coverageMerged_{ED}.union.bed"),
        partition=temp(cov_dir + "mosdepth_by_threshold/coverageMerged_{ED}.partition.bed"),
        combined=temp(cov_dir + "mosdepth_by_threshold/coverageMerged_{ED}.combined.bed"),
        out=cov_dir + "gencov.mismatch.{ED}.out"
    conda: 
        "../envs/bedops.yaml"
    shell:
        """
        bedops --everything {input.bed_hetero} {input.bed_homo} > {output.union}
        bedops --partition {output.union} > {output.partition}
        bedmap --echo --echo-map-id --delim '\t' {output.partition} {output.union} | sed 's/;/\t/g' | awk '$3-$2=="5000" {{print}}' > {output.combined}
        bedops --element-of 1 {output.combined} {input.windows} > {output.out}
        """


rule normalize_cov_mean:
    input:
        cov_dir + "mosdepth_by_threshold/coverageMerged_mismatch.{ED}.bed"
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
        python3 workflow/scripts/normalize_genCov.py {input} no-synteny {params.hetero} {params.homo} > {output}
        """
