# Written by Hanna Sigeman, Nov 2019

##########################################################
##################### INDEX GENOME #######################
##########################################################

rule index_fasta_bwa:
    input:
        REF_FASTA
    output:
        REF_FASTA + ".bwt"
    priority : 80
    message:
        """--- Indexing {input} with BWA index."""
    threads: 2
    shell:
        """
        bwa index {input}
        """

rule index_fasta_samtools:
    input:
        REF_FASTA
    output:
        REF_FASTA + ".fai"
    priority : 70
    threads: 2
    shell:
        """
        samtools faidx {input}
        """

##########################################################
######################## MAPPING #########################
##########################################################

rule map:
    input:
        R1= FQ_DIR + "{S}_forward_paired.fq.gz",
        R2= FQ_DIR + "{S}_reverse_paired.fq.gz",
        ref = REF_FASTA,
        ref_bwt = REF_FASTA + ".bwt"
    output:
        temp(MAP_DIR + "{S}.bam")
    message: "Mapping reads to ref"
    threads: 20
    params:
        rg="\"@RG\\tID:{S}\\tSM:{S}\""
    shell:
        """
        bwa mem -t {threads} -M -R {params.rg} {input.ref} {input.R1} {input.R2} | samtools view -Sb - > {output}
        """

rule sort_bam:
    input:
        MAP_DIR + "{S}.bam"
    output:
        temp(MAP_DIR + "{S}.sorted.bam")
    log: MAP_DIR + "{S}.sorted.status"
    threads: 3
    params:
        tmpdir = MAP_DIR + "{S}_temp_sort/"
    shell:
        """
        mkdir {params.tmpdir}
        samtools sort -@ {threads} {input} -T {params.tmpdir} > {output}
        rm -r {params.tmpdir}
        echo "DONE" > {log}
        """

rule remove_duplicates:
    input:
        MAP_DIR + "{S}.sorted.bam"
    output:
        MAP_DIR + "{S}.sorted.nodup.nm.all.bam"
    log: MAP_DIR + "{S}.sorted.nodup.nm.all.status"
    params:
        tmpdir = MAP_DIR + "{S}_temp_dupl/"
    shell:
        """
        mkdir {params.tmpdir}
        picard MarkDuplicates -Xmx10g MAX_FILE_HANDLES=500 REMOVE_DUPLICATES=true I={input} O={output} M={input}_duplicatedata.txt TMP_DIR={params.tmpdir}
        rm -r {params.tmpdir}
        echo "DONE" > {log}
        """

rule index_bam:
    input:
        MAP_DIR + "{S}.sorted.nodup.nm.{ED}.bam"
    output:
        MAP_DIR + "{S}.sorted.nodup.nm.{ED}.bam.bai"
    threads: 1
    shell:
        """
        samtools index {input}
        """

rule mismatch_bam:
    input:
        bam = MAP_DIR + "{S}.sorted.nodup.nm.all.bam",
        ref = REF_FASTA,
        bai = MAP_DIR + "{S}.sorted.nodup.nm.all.bam.bai"
    output:
        MAP_DIR + "{S}.sorted.nodup.nm.0.{ED, [0-9]+}.bam"
    threads: 2
    params:
        "\"NM:i:[0-{ED}]\""
    shell:
        """
        samtools view -@ {threads} {input.bam} | grep -E {params} | samtools view -bS -T {input.ref} - > {output}
        """

rule samtools_stats:
    input:
        MAP_DIR + "{S}.sorted.nodup{ED}bam"
    output:
        MAP_DIR + "{S}.sorted.nodup{ED}stats"
    shell:
        """
        samtools stats {input} > {output}
        """

##########################################################
#################### GENOME COVERAGE #####################
##########################################################

rule gencov_prepare_fasta:
    input:
        REF_FASTA + ".fai"
    output:
        filter_fai = GENCOV_DIR_REF + REF_SPECIES + ".filter." + MIN_SIZE_SCAFFOLD + ".fasta.fai",
        windows = GENCOV_DIR_REF + "genome_5kb_windows.out"
    params:
        MIN_SIZE_SCAFFOLD
    threads: 1
    shell:
        """
	cat {input} | awk '$2>= {params} {{print $0}}' > {output.filter_fai}
        bedtools makewindows -g {output.filter_fai} -w 5000 -s 5000 | awk '$3-$2==5000 {{print}}' > {output.windows}
        """

rule gencov_bedtools:
    input:
        bam_hetero = expand(MAP_DIR + "{heterogametic}.sorted.nodup.nm.{{ED}}.bam", heterogametic = HETEROGAMETIC),
        bai_hetero = expand(MAP_DIR + "{heterogametic}.sorted.nodup.nm.{{ED}}.bam.bai", heterogametic = HETEROGAMETIC),
        bam_homo = expand(MAP_DIR + "{homogametic}.sorted.nodup.nm.{{ED}}.bam", homogametic = HOMOGAMETIC),
        bai_homo = expand(MAP_DIR + "{homogametic}.sorted.nodup.nm.{{ED}}.bam.bai", homogametic = HOMOGAMETIC),
        bed = GENCOV_DIR_REF + "genome_5kb_windows.out"
    output:
        protected(GENCOV_DIR + SPECIES + ".gencov.nodup.nm.{ED}.out")
    threads: 2
    shell:
        """
        bedtools multicov -bams {input.bam_hetero} {input.bam_homo} -bed {input.bed} -p -q 20 > {output}
        """

##########################################################
#################### VARIANT CALLING #####################
##########################################################

rule freebayes_prep:
    input:
        REF_FASTA + ".fai"
    output:
        filter_fai = VCF_DIR_REF + REF_SPECIES + ".filter." + MIN_SIZE_SCAFFOLD + ".fasta.fai",
        regions = VCF_DIR_REF + REF_SPECIES + ".100kbp.regions"
    params:
        MIN_SIZE_SCAFFOLD
    threads: 4
    shell:
        """
	cat {input} | awk '$2>= {params} {{print $0}}' > {output.filter_fai}
        fasta_generate_regions.py {output.filter_fai} 100000 > {output.regions}
        """


rule freebayes_parallel:
    input:
        ref = REF_FASTA,
        regions = VCF_DIR_REF + REF_SPECIES + ".100kbp.regions",
        samples = expand(MAP_DIR + "{S}.sorted.nodup.nm.all.bam", S = ID),
	    bai = expand(MAP_DIR + "{S}.sorted.nodup.nm.all.bam.bai", S = ID)
    output:
        vcf = temp(VCF_DIR + SPECIES + ".vcf"),
        gz = VCF_DIR + SPECIES + ".vcf.gz"
    log: VCF_DIR + SPECIES + ".vcf.log"
    priority : 60
    threads: 18
    params:
        tmpdir=VCF_DIR + "temp/"
    shell:
        """
        mkdir {params.tmpdir}
        export TMPDIR={params.tmpdir}
        freebayes-parallel {input.regions} {threads} -f {input.ref} {input.samples} > {output.vcf}
        rm -r {params.tmpdir}

        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}

        echo "DONE" > {log}
        """

rule vcftools_filter:
    input:
        VCF_DIR + SPECIES + ".vcf.gz"
    output:
        vcf = temp(VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf"),
        gz = VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf.gz"
    shell:
        """
        vcftools --gzvcf {input} --min-alleles 2 --max-alleles 2 --remove-filtered-geno-all --minQ 20 --minDP 3 --recode --stdout > {output.vcf}

        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        """

##########################################################
##################### SNP ANALYSIS #######################
##########################################################

rule proportion_heterozygosity:
    input:
        VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf.gz"
    output:
        diff_het = temp(VCF_DIR + SPECIES + ".diffHeterozygosity.bed"),
        het = temp(VCF_DIR + SPECIES + ".heterozygosity.bed"),
    log: VCF_DIR + SPECIES + ".proportion_heterozygosity.log"
    params:
        hetero = expand("het:{heterogametic}", heterogametic = HETEROGAMETIC),
        homo = expand("homo:{homogametic}", homogametic = HOMOGAMETIC)
    threads: 1
    shell:
        """
        python3 code/calculate_hetDiff.py {input} {output.diff_het} {params.hetero} {params.homo} > {log}

        python3 code/heterozygosity_per_indv.py {input} {output.het} {params.hetero} {params.homo} >> {log}
        """

rule proportion_heterozygosity_window:
    input:
        diff_het_sorted = VCF_DIR + SPECIES + ".diffHeterozygosity.bed",
        het_sorted = VCF_DIR + SPECIES + ".heterozygosity.bed",
        windows = GENCOV_DIR_REF + "genome_5kb_windows.out"
    output:
        diff_het_sorted_window = temp(VCF_DIR + SPECIES + ".diffHeterozygosity.5kb.windows.bed"),
        diff_het_sorted_window_mean = VCF_DIR + SPECIES + ".diffHeterozygosity.5kb.windows.mean.bed",
        het_sorted_window = temp(VCF_DIR + SPECIES + ".heterozygosity.5kb.windows.bed"),
        het_sorted_window_mean = VCF_DIR + SPECIES + ".heterozygosity.5kb.windows.mean.bed"
    threads: 1
    shell:
        """
        bedtools intersect -a {input.windows} -b {input.diff_het_sorted} -wa -wb | cut -f 1-3,7 > {output.diff_het_sorted_window}
        cat {output.diff_het_sorted_window} | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' | awk '{{count[$1]++; sum[$1]+=$2}} END {{for(i in count) {{m = sum[i]/count[i]; print i, m}}}}' |  sed 's/STARTCOORD/\t/' | sed 's/ENDCOORD/\t/' | sed 's/ /\t/g' > {output.diff_het_sorted_window_mean}

        bedtools intersect -a {input.windows} -b {input.het_sorted} -wa -wb | cut -f 1-3,7- > {output.het_sorted_window}
        cat {output.het_sorted_window} | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' | awk ' {{c[$1]++; for (i=2;i<=NF;i++) {{ s[$1"."i]+=$i}}; }} END {{for (k in c) {{printf "%s\t", k; for(i=2;i<NF;i++) printf "%.1f\\t", s[k"."i]/c[k]; printf "%.1f\\n", s[k"."NF]/c[k];}}}}' |  sed 's/STARTCOORD/\t/' | sed 's/ENDCOORD/\t/' | sed 's/ /\t/g' > {output.het_sorted_window_mean}
	"""

#rule allele_frequency:
#    input:
#        VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf.gz"
#    output:
#        hetero = VCF_DIR + SPECIES + ".allFreq.heterogametic.out",
#        homo = VCF_DIR + SPECIES + ".allFreq.homogametic.out"
#    params:
#        hetero = expand("--indv {heterogametic}", heterogametic = HETEROGAMETIC),
#        homo = expand("--indv {homogametic}", homogametic = HOMOGAMETIC)
#    shell:
#        """
#        vcftools --gzvcf {input} {params.hetero} --freq --stdout > {output.hetero}
#        vcftools --gzvcf {input} {params.homo} --freq --stdout > {output.homo}
#        """

#rule filter_allele_frequency:
#    input:
#        hetero = VCF_DIR + SPECIES + ".allFreq.heterogametic.out",
#        homo = VCF_DIR + SPECIES + ".allFreq.homogametic.out"
#    output:
#        bed = temp(VCF_DIR + SPECIES + ".allFreq.bed"),
#        sorted = VCF_DIR + SPECIES + ".allFreq.sorted.bed"
#    threads: 1
#    shell:
#        """
#        python3 code/filter_allFreq.py {input.hetero} heterogametic > {output.bed}
#        python3 code/filter_allFreq.py {input.homo} homogametic >> {output.bed}
#        sort {output.bed} -k 1,1 -k 2,2 > {output.sorted}
#        """

##########################################################
####################### RESULTS ##########################
##########################################################

rule plotting:
    input:
        cov = expand(RESULTDIR + SPECIES + "{{synteny}}{{gencov}}.nodup.nm.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = RESULTDIR + SPECIES + "{synteny}diffHeterozygosity.{bp}bp.out"
    output:
        touch(RESULTDIR + SPECIES + "{synteny}{gencov}.plotting.{bp}bp.done")
    threads: 1
    params:
        out_circlize = protected(RESULTDIR + SPECIES + "{synteny}{gencov}.circlize.{bp}bp.pdf"),
        out_scatter = protected(RESULTDIR + SPECIES + "{synteny}{gencov}.scatter.{bp}bp.pdf"),
        chromosomes = CHROMOSOMES
    wildcard_constraints:
        gencov = "gencov|N1"
    shell:
        """
        Rscript code/plot_windows.R {input.cov} {input.snp} {params.out_circlize} {params.out_scatter} {params.chromosomes}
        """

rule plotting_chr:
    input:
        cov = expand(RESULTDIR + SPECIES + "{{synteny}}gencov.nodup.nm.{ED}.chr.out", ED = EDIT_DIST),
        snp = RESULTDIR + SPECIES + "{synteny}diffHeterozygosity.chr.out"
    output:
        touch(RESULTDIR + SPECIES + "{synteny}plotting_chr.done")
    params:
        out_scatter2D = protected(RESULTDIR + SPECIES + "{synteny}chr_scatter2D.pdf"),
        out_scatter3D = protected(RESULTDIR + SPECIES + "{synteny}chr_scatter3D.pdf"),
        chromosomes = CHROMOSOMES
    threads: 1
    shell:
        """
        Rscript code/scatterplot_chr.R {input.cov} {input.snp} {params.out_scatter2D} {params.out_scatter3D} {params.chromosomes}
        """

##########################################################
################### MODIFY REF GENOME ####################
##########################################################

#rule modify_genome:
#    input:
#        vcf = VCF_DIR + SPECIES + ".vcf.gz",
#        ref = REF_FASTA
#    output:
#        noN_vcf = temp(VCF_DIR + SPECIES + ".noN.vcf.gz"),
#        vcf = temp(VCF_DIR + SPECIES + ".non-ref-af_05_biallelic_qual.vcf"),
#        gz = VCF_DIR + SPECIES + ".non-ref-af_05_biallelic_qual.vcf.gz",
#        ref = REF_DIR + REF_NAME + "_nonRefAf_consensus.fasta"
#    shell:
#        """
#        zcat {input.vcf} | awk '$4 !~ /N/' | gzip -c > {output.noN_vcf}
#        vcftools --gzvcf {output.noN_vcf} --non-ref-af 0.5 --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout > {output.vcf}
#        bgzip -c {output.vcf} > {output.gz}
#        tabix -p vcf {output.gz}
#        cat {input.ref} | bcftools consensus {output.gz} > {output.ref}
#        """

rule modify_genome:
    input:
        vcf = VCF_DIR + SPECIES + ".vcf.gz",
        ref = REF_FASTA
    output:
        vcf = temp(VCF_DIR + SPECIES + ".non-ref-af_05_biallelic_qual.vcf"),
        gz = VCF_DIR + SPECIES + ".non-ref-af_05_biallelic_qual.vcf.gz",
        ref = REF_DIR + REF_NAME + "_nonRefAf_consensus.fasta"
    shell:
        """
        vcftools --gzvcf {input.vcf} --non-ref-af 0.5 --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        cat {input.ref} | bcftools consensus {output.gz} > {output.ref}
        """
