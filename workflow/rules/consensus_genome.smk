rule modify_genome:
    input:
        vcf = outdir + "variant_calling/" + ref_genome_name_simple + ".vcf.gz",
        ref = ref_genome
    output:
        vcf = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".non-ref-af_05_biallelic_qual.vcf"),
        gz = outdir + "variant_calling/" + ref_genome_name_simple + ".non-ref-af_05_biallelic_qual.vcf.gz",
        ref = outdir + "consensus_genome/" + ref_genome_name_simple + "_nonRefAf_consensus.fasta"
    shell:
        """
        vcftools --gzvcf {input.vcf} --non-ref-af 0.5 --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        cat {input.ref} | bcftools consensus {output.gz} > {output.ref}
        """
