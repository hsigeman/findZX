rule freebayes_prep:
    input:
        fai = ref_genome + ".fai",
        samples = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.all.bam", u=units.itertuples()),
	samples_bai = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.all.bam.bai", u=units.itertuples()),
    output:
        filter_fai = outdir + "variant_calling/" + ref_genome_name_simple + ".filter." + MIN_SIZE_SCAFFOLD + ".list",
        regions = outdir + "variant_calling/" + ref_genome_name_simple + ".freebayes.regions",
        regions_filter = outdir + "variant_calling/" + ref_genome_name_simple + ".freebayes.regions.filter"
    params:
        MIN_SIZE_SCAFFOLD
    conda: 
        "../envs/python_gawk.yaml"
    shell:
        """
	    cat {input.fai} | awk '$2>= {params} {{print $1}}' > {output.filter_fai}
        python3 code/split_ref_by_bai_datasize.py {input.samples} -r {input.fai} > {output.regions}
        join {output.regions} {output.filter_fai} | sed 's/ /\t/g' | sed 's/\t/:/' | sed 's/\t/-/' > {output.regions_filter}
        """


rule freebayes:
    input:
        ref = ref_genome,
        samples = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.all.bam", u=units.itertuples()),
        indexes=expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.all.bam.bai", u=units.itertuples()),
        regions_filter=outdir + "variant_calling/" + ref_genome_name_simple + ".freebayes.regions.filter"
    output:
        vcf = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".vcf"),
        log = outdir + "variant_calling/freebayes_done.log"
    params:
        extra="--use-best-n-alleles 4 --skip-coverage 1000",         # optional parameters
    log:
        logs_dir + "freebayes/freebayes.log"
    threads: threads_max
    conda: 
        "../envs/freebayes.yaml"
    shell:
        """
        freebayes-parallel {input.regions_filter} {threads} {params.extra} -f {input.ref} {input.samples} > {output.vcf}
        echo "DONE" > {output.log}
        """


rule bgzip_tabix:
    input:
        vcf = outdir + "variant_calling/" + ref_genome_name_simple + ".vcf",
        log = outdir + "variant_calling/freebayes_done.log"
    output:
        outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.gz"
    log: outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.log"
    conda: 
        "../envs/vcftools_filter.yaml"
    shell:
        """
        bgzip -c {input.vcf} > {output}
        tabix -p vcf {output}
        """


rule vcftools_filter:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.gz"
    output:
        vcf = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf"),
        gz = outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    conda: 
        "../envs/vcftools_filter.yaml"
    shell:
        """
        vcftools --gzvcf {input} --min-alleles 2 --max-alleles 2 --remove-filtered-geno-all --minQ 20 --minDP 3 --recode --stdout > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        """