rule platypus:
    input:
        ref = ref_genome,
        samples = expand(dedup_dir + "{u.sample}__{u.group}.sorted.dedup.mismatch.unfiltered.bam", u=units.itertuples()),
        indexes=expand(dedup_dir + "{u.sample}__{u.group}.sorted.dedup.mismatch.unfiltered.bam.bai", u=units.itertuples()),
    output:
        vcf = temp(vcf_dir + ref_genome_name_simple + ".vcf"),
        log = vcf_dir + "freebayes_done.log"
    params:
        files=lambda wildcards, input: ','.join(input.samples)
    log:
        logs_dir + "platypus/v.log"
    threads: threads_max
    conda: 
        "../envs/platypus.yaml"
    message:
        "Variant calling (platypus-variant)"
    shell:
        """
        platypus callVariants --bamFiles={params.files} --refFile={input.ref} --output={output.vcf} --nCPU={threads}
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

rule het_calc:
    input:
        vcf_dir + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    output:
        vcf_dir + ref_genome_name_simple + ".biallelic.minQ20.minDP3.het",
    conda: 
        "../envs/vcftools_filter.yaml"
    message:
        "Filter VCF file"
    shell:
        """
        vcftools --gzvcf {input} --het --stdout > {output}
        """

rule het_calc_genome:
    input:
        het=vcf_dir + ref_genome_name_simple + ".biallelic.minQ20.minDP3.het",
        assembly_stats=qc_dir + "assembly_stats/" + ref_genome_name_simple + "_stats.txt",
    output:
        report(vcf_dir + ref_genome_name_simple + ".heterozygosity.perc.csv", category="4. Stats"),
    message:
        "Calculate heterozygosity"
    shell:
        """
        length=$(cat {input.assembly_stats} | grep total_length | cut -f 3) 
        echo "Sample,Heterozygous_sites,Genome_length,Percentage_heterozygosity" > {output}
        cat {input.het} | grep -v INDV | awk -v var="$length" '{{print $1,($4-$2),var,($4-$2)/var*100}}' | sed 's/ /,/g' >> {output}
        """
