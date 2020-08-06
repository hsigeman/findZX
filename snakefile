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
    log: MAP_DIR + "{S}.sorted.log"
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
        protected(MAP_DIR + "{S}.sorted.nodup.nm.all.bam")
    log: MAP_DIR + "{S}.sorted.nodup.nm.all.log"
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
        ref = REF_FASTA
    output: 
        MAP_DIR + "{S}.sorted.nodup.nm.0.{ED, [0-9]+}.bam"
    threads: 2
    params:
        "\"NM:i:[0-{ED}]\""
    shell: 
        """
        samtools view -@ {threads} {input.bam} | grep -E {params} | samtools view -bS -T {input.ref} - > {output}
        """

rule flagstat:
    input:
        MAP_DIR + "{S}.sorted{V}bam"
    output:
        MAP_DIR + "{S}.sorted{V}flagstat"
    shell:
        """
        samtools flagstat {input} > {output}  
        """

##########################################################  
#################### GENOME COVERAGE #####################       
########################################################## 

rule gencov_prepare_fasta:
    input: 
        REF_FASTA + ".fai"
    output: 
        GENCOV_DIR_REF + "genome_5kb_windows.out"
    threads: 1
    shell: 
        """
        bedtools makewindows -g {input} -w 5000 -s 5000 > {output}
        """

rule gencov_bedtools:
    input: 
        bam_hetero = expand(MAP_DIR + "{heterogametic}.sorted.nodup.nm.{{V}}.bam", heterogametic = HETEROGAMETIC),
        bai_hetero = expand(MAP_DIR + "{heterogametic}.sorted.nodup.nm.{{V}}.bam.bai", heterogametic = HETEROGAMETIC),
        bam_homo = expand(MAP_DIR + "{homogametic}.sorted.nodup.nm.{{V}}.bam", homogametic = HOMOGAMETIC),
        bai_homo = expand(MAP_DIR + "{homogametic}.sorted.nodup.nm.{{V}}.bam.bai", homogametic = HOMOGAMETIC),
        bed = GENCOV_DIR_REF + "genome_5kb_windows.out"
    output: 
        protected(GENCOV_DIR + SPECIES + ".gencov.nodup.nm.{V}.out")
    threads: 2
    shell: 
        """
        bedtools multicov -bams {input.bam_hetero} {input.bam_homo} -bed {input.bed} > {output}
        """

##########################################################  
#################### VARIANT CALLING #####################      
########################################################## 

rule freebayes_prep:
    input: 
        REF_FASTA + ".fai"
    output: 
        VCF_DIR_REF + REF_SPECIES + ".100kbp.regions"
    threads: 4
    shell: 
        """
        fasta_generate_regions.py {input} 100000 > {output}
        """

rule freebayes_parallel:
    input: 
        ref = REF_FASTA,
        regions = VCF_DIR_REF + REF_SPECIES + ".100kbp.regions",
        samples = expand(MAP_DIR + "{S}.sorted.nodup.nm.all.bam", S = ID)
    output: 
        vcf = temp(VCF_DIR + SPECIES + ".vcf"),
        gz = protected(VCF_DIR + SPECIES + ".vcf.gz")
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
        gz = protected(VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf.gz")
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
        diff_het_sorted = VCF_DIR + SPECIES + ".diffHeterozygosity.sorted.bed",
        het_sorted = VCF_DIR + SPECIES + ".heterozygosity.sorted.bed"
    log: VCF_DIR + SPECIES + ".proportion_heterozygosity.log"
    params:
        hetero = expand("het:{heterogametic}", heterogametic = HETEROGAMETIC),
        homo = expand("homo:{homogametic}", homogametic = HOMOGAMETIC)
    threads: 1
    shell:
        """
        python3 code/calculate_hetDiff.py <(less {input}) {output.diff_het} {params.hetero} {params.homo}
        sort {output.diff_het} -k 1,1 -k 2,2 > {output.diff_het_sorted} > {log}

        python3 code/heterozygosity_per_indv.py <(less {input}) {output.het} {params.hetero} {params.homo}
        sort {output.het} -k 1,1 -k 2,2 > {output.het_sorted} > {log}
        """

rule allele_frequency:
    input:
        VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf.gz"
    output:
        hetero = VCF_DIR + SPECIES + ".allFreq.heterogametic.out",
        homo = VCF_DIR + SPECIES + ".allFreq.homogametic.out"
    params:
        hetero = expand("--indv {heterogametic}", heterogametic = HETEROGAMETIC),
        homo = expand("--indv {homogametic}", homogametic = HOMOGAMETIC)
    shell:
        """
        vcftools --gzvcf {input} {params.hetero} --freq --stdout > {output.hetero}
        vcftools --gzvcf {input} {params.homo} --freq --stdout > {output.homo}
        """

rule filter_allele_frequency:
    input:
        hetero = VCF_DIR + SPECIES + ".allFreq.heterogametic.out",
        homo = VCF_DIR + SPECIES + ".allFreq.homogametic.out"
    output:
        bed = temp(VCF_DIR + SPECIES + ".allFreq.bed"),
        sorted = VCF_DIR + SPECIES + ".allFreq.sorted.bed"
    threads: 1
    shell:
        """
        python3 code/filter_allFreq.py {input.hetero} heterogametic > {output.bed}
        python3 code/filter_allFreq.py {input.homo} homogametic >> {output.bed}
        sort {output.bed} -k 1,1 -k 2,2 > {output.sorted}
        """

##########################################################
####################### RESULTS ##########################
##########################################################

rule print_report:
    input:
        cov0 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.0.chr.out",
        cov2 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.2.chr.out",
        cov4 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.4.chr.out",
        snp = RESULTDIR + SPECIES + "{synteny}diffHeterozygosity.chr.out"
    output:
        stats = RESULTDIR + SPECIES + "{synteny}statistics.tsv",
        html = RESULTDIR + SPECIES + "{synteny}report.html"
    params:
        species = SPECIES,
        ref = REF_SPECIES,
        heterogametic = expand("het:{het}", het = HETEROGAMETIC),
        homogametic = expand("homo:{hom}", hom =  HOMOGAMETIC),
        synt = "{synteny}"
    threads: 1
    shell:
        """
        join -1 1 -2 1 {input.cov0} {input.cov2} | join -1 1 -2 1 - {input.cov4} | -1 1 -2 1 - {input.snp} | sed 's/ /\t/g' | cut -f 1,2,4,6,8,9 > {output.stats}

        python3 code/make_info_html.py {params.species} {params.ref} {output.stats} {params.synt} {params.heterogametic} {params.homogametic} > {output.html}
        """

rule plotting:
    input: 
        cov0 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.0.1Mbp.out",
        cov2 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.2.1Mbp.out",
        cov4 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.4.1Mbp.out",
        snp = RESULTDIR + SPECIES + "{synteny}diffHeterozygosity.1Mbp.out"
    output: 
        circlize = protected(RESULTDIR + SPECIES + "{synteny}circlize.pdf"),
        scatter = protected(RESULTDIR + SPECIES + "{synteny}scatter.pdf")
    threads: 1
    params:
        chromosomes = CHOMOSOMES
    shell: 
        """
        Rscript code/plot_windows.R {input.cov0} {input.cov2} {input.cov4} {input.snp} {output.circlize} {output.scatter} {params.chromosomes}
        """

rule plotting_chr:
    input:
        cov0 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.0.chr.out",
        cov2 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.2.chr.out",
        cov4 = RESULTDIR + SPECIES + "{synteny}gencov.nodup.nm.0.4.chr.out",
        snp = RESULTDIR + SPECIES + "{synteny}diffHeterozygosity.chr.out"
    output:
        scatter2D = protected(RESULTDIR + SPECIES + "{synteny}chr_scatter2D.pdf"),
        scatter3D = protected(RESULTDIR + SPECIES + "{synteny}chr_scatter3D.pdf")
    threads: 1
    shell:
        """
        Rscript code/scatterplot_chr.R {input.cov0} {input.cov2} {input.cov4} {input.snp} {output.scatter2D} {output.scatter3D}
        """

##########################################################
################### MODIFY REF GENOME ####################
##########################################################

rule modify_genome:
    input:
        vcf = VCF_DIR + SPECIES + ".vcf.gz",
        ref = REF_FASTA
    output:
        vcf = temp(VCF_DIR + SPECIES + ".non-ref-ac_05_biallelic_qual.vcf"),
        gz = VCF_DIR + SPECIES + ".non-ref-ac_05_biallelic_qual.vcf.gz",
        ref = REF_DIR + REF_NAME + "_nonRefAf_consensus.fasta"
    shell:
        """
        vcftools --vcf {input.vcf} --non-ref-af 0.5 --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        cat {input.ref} | bcftools consensus -I {output.gz} > {output.ref}
        """

