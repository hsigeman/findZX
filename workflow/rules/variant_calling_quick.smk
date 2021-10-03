rule add_ID_to_VCF:
    input:
        vcf_dir + "{sample}__{group}.vcf"
    output:
        vcf_dir + "{sample}__{group}.ID.vcf"
    params:
        get_ID
    shell:
        """
        sed 's/unknown/{params}/' {input} | grep "^#" | sed -e "s/Number=[A-Z]/Number=1/" > {output}
        grep -v "#" {input} | grep PASS | sed -e "s/RC=[0-9];//" | sed -e "s/RC=[0-9]+;//" >> {output}
        """

rule bgzip_tabix:
    input:
        vcf_dir + "{sample}__{group}.ID.vcf"
    output:
        vcf_dir + "{sample}__{group}.vcf.gz"
#   log: vcf_dir + ref_genome_name_simple + ".vcf.log"
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
        bam_hetero = expand(vcf_dir + "{u.sample}__{u.group}.vcf.gz", zip, u=heterogametic.itertuples()),
        bam_homo = expand(vcf_dir + "{u.sample}__{u.group}.vcf.gz", zip, u=homogametic.itertuples()),
    output:
        vcf= vcf_dir + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    conda: 
        "../envs/vcftools_filter.yaml"
    shell: 
        "vcf-merge {input.bam_hetero} {input.bam_homo} -R \"0/0\" | bgzip -c > {output}"
