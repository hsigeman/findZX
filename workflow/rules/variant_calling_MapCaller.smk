rule add_ID_to_VCF:
    input:
        outdir + "variant_calling/" + "{sample}__{group}.vcf"
    output:
        outdir + "variant_calling/" + "{sample}__{group}.ID.vcf"
    params:
        get_ID
    shell:
        """
        sed 's/unknown/{params}/' {input} | grep "^#" | sed -e "s/Number=[A-Z]/Number=1/" > {output}
        grep -v "#" {input} | grep PASS | sed -e "s/RC=[0-9];//" | sed -e "s/RC=[0-9]+;//" >> {output}
        """

rule bgzip_tabix:
    input:
        outdir + "variant_calling/" + "{sample}__{group}.ID.vcf"
    output:
        outdir + "variant_calling/" + "{sample}__{group}.vcf.gz"
 #   log: outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.log"
    conda: 
        "../envs/vcftools_filter.yaml"
    message:
        "Compress VCF file"
    shell:
        """
        bgzip -c {input} > {output}
        tabix -p vcf {output}
        """

rule merge_vcf:
    input: 
        bam_hetero = expand(outdir + "variant_calling/" + "{u.sample}__{u.group}.vcf.gz", zip, u=heterogametic.itertuples()),
        bam_homo = expand(outdir + "variant_calling/" + "{u.sample}__{u.group}.vcf.gz", zip, u=homogametic.itertuples()),
    output:
        vcf=outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.gz"
    conda: 
        "../envs/modify_genome.yaml"
    shell: 
        "bcftools merge {input.bam_hetero} {input.bam_homo} -Oz -o {output}"


rule vcftools_filter:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.gz"
    output:
        vcf = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf"),
        gz = outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    conda: 
        "../envs/vcftools_filter.yaml"
    message:
        "Filter VCF file"
    shell:
        """
        vcftools --gzvcf {input} --recode --stdout > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        """