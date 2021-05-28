rule freebayes_prep:
    input:
        fai = ref_genome + ".fai",
        samples = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.all.bam", u=units.itertuples()),
	samples_bai = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.all.bam.bai", u=units.itertuples()),
    output:
        filter_fai = outdir + "variant_calling/" + ref_genome_name_simple + ".filter." + MIN_SIZE_SCAFFOLD + ".fasta.fai",
        regions = outdir + "variant_calling/" + ref_genome_name_simple + ".freebayes.regions"
    params:
        MIN_SIZE_SCAFFOLD
    shell:
        """
	    cat {input.fai} | awk '$2>= {params} {{print $0}}' > {output.filter_fai}
        python3 code/split_ref_by_bai_datasize.py {input.samples} -r {input.fai} | sed 's/ /\t/g' | sed 's/\t/:/' | sed 's/\t/-/'  > {output.regions}
        """


rule freebayes:
    input:
        ref = ref_genome,
        samples = expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.all.bam", u=units.itertuples()),
        indexes=expand(dedup_dir + "{u.sample}__{u.unit}.sorted.dedup.nm.all.bam.bai", u=units.itertuples()),
        regions=outdir + "variant_calling/" + ref_genome_name_simple + ".freebayes.regions"
    output:
        temp(outdir + "variant_calling/" + ref_genome_name_simple + ".vcf")
    params:
        extra="--use-best-n-alleles 4 --skip-coverage 1000",         # optional parameters
    log:
        logs_dir + "freebayes/freebayes.log"
    threads: threads_max
    shell:
        """
        freebayes-parallel {input.regions} {threads} {params.extra} -f {input.ref} {input.samples} > {output}
        """


rule bgzip_tabix:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".vcf"
    output:
        outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.gz"
    log: outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.log"
    shell:
        """
        bgzip -c {input} > {output}
        tabix -p vcf {output}
        """


rule vcftools_filter:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.gz"
    output:
        vcf = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf"),
        gz = outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    shell:
        """
        vcftools --gzvcf {input} --min-alleles 2 --max-alleles 2 --remove-filtered-geno-all --minQ 20 --minDP 3 --recode --stdout > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        """