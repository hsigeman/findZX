rule freebayes_prep:
    input:
        fai = ref_genome + ".fai",
        samples = expand(dedup_dir + "{u.sample}__{u.group}.sorted.dedup.mismatch.unfiltered.bam", u=units.itertuples()),
	samples_bai = expand(dedup_dir + "{u.sample}__{u.group}.sorted.dedup.mismatch.unfiltered.bam.bai", u=units.itertuples()),
    output:
        filter_fai = vcf_dir + ref_genome_name_simple + ".filter." + MIN_SIZE_SCAFFOLD + ".list",
        regions = vcf_dir + ref_genome_name_simple + ".freebayes.regions",
        regions_filter = vcf_dir + ref_genome_name_simple + ".freebayes.regions.filter"
    params:
        MIN_SIZE_SCAFFOLD
    log:
        logs_dir + "freebayes/freebayes_prep.log"
    conda: 
        "../envs/python_gawk.yaml"
    message:
        "Prepare region file for variant calling (freebayes)"
    shell:
        """
	    cat {input.fai} | awk '$2>= {params} {{print $1}}' | sort -k1,1 > {output.filter_fai}
        python3 workflow/scripts/split_ref_by_bai_datasize.py {input.samples} -r {input.fai} | sed 's/ /\t/g' | bedtools sort > {output.regions} 2> {log}
        join <(sort {output.regions}) <(sort {output.filter_fai}) | sed 's/ /\t/g' | sed 's/\t/:/' | sed 's/\t/-/' > {output.regions_filter}
        """


rule freebayes:
    input:
        ref = ref_genome,
        samples = expand(dedup_dir + "{u.sample}__{u.group}.sorted.dedup.mismatch.unfiltered.bam", u=units.itertuples()),
        indexes=expand(dedup_dir + "{u.sample}__{u.group}.sorted.dedup.mismatch.unfiltered.bam.bai", u=units.itertuples()),
        regions_filter= vcf_dir + ref_genome_name_simple + ".freebayes.regions.filter"
    output:
        vcf = temp(vcf_dir + ref_genome_name_simple + ".vcf"),
        log = vcf_dir + "freebayes_done.log"
    params:
        extra="--use-best-n-alleles 4 --skip-coverage 1000",         # optional parameters
    log:
        logs_dir + "freebayes/freebayes.log"
    threads: threads_max
    conda: 
        "../envs/freebayes.yaml"
    message:
        "Variant calling (freebayes-parallel)"
    shell:
        """
        freebayes-parallel {input.regions_filter} {threads} {params.extra} -f {input.ref} {input.samples} > {output.vcf}
        echo "DONE" > {output.log}
        """


rule bgzip_tabix:
    input:
        vcf = vcf_dir + ref_genome_name_simple + ".vcf",
        log = vcf_dir + "freebayes_done.log"
    output:
        vcf_dir + ref_genome_name_simple + ".vcf.gz"
    log: vcf_dir + ref_genome_name_simple + ".vcf.log"
    conda: 
        "../envs/vcftools_filter.yaml"
    message:
        "Compress VCF file"
    shell:
        """
        bgzip -c {input.vcf} > {output}
        tabix -p vcf {output}
        """


rule vcftools_filter:
    input:
        vcf_dir + ref_genome_name_simple + ".vcf.gz"
    output:
        vcf = temp(vcf_dir + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf"),
        gz = vcf_dir + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    conda: 
        "../envs/vcftools_filter.yaml"
    message:
        "Filter VCF file"
    shell:
        """
        vcftools --gzvcf {input} --min-alleles 2 --max-alleles 2 --remove-filtered-geno-all --minQ 20 --minDP 3 --recode --stdout > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        """