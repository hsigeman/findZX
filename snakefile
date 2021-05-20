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
        mkdir -p {params.tmpdir}
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
        filter_fai = GENCOV_DIR_REF + REF_NAME + ".filter." + MIN_SIZE_SCAFFOLD + ".fasta.fai",
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
        GENCOV_DIR + SPECIES + ".gencov.nodup.nm.{ED}.out"
    threads: 2
    shell:
        """
        bedtools multicov -bams {input.bam_hetero} {input.bam_homo} -bed {input.bed} -p -q 20 > {output}
        """
	
rule normalize_cov_mean:
    input:
        GENCOV_DIR + SPECIES + ".gencov.nodup.nm.{ED}.out"
    output:
        GENCOV_DIR + SPECIES + ".gencov.nodup.nm.{ED}.norm.sexAverage.out"
    threads: 1
    params:
        hetero = expand("het:{heterogametic}", heterogametic = HETEROGAMETIC),
        homo = expand("homo:{homogametic}", homogametic = HOMOGAMETIC)
    shell:
        """
        python3 code/normalize_genCov.py {input} no-synteny {params.hetero} {params.homo} > {output}
        """


##########################################################
#################### VARIANT CALLING #####################
##########################################################

rule freebayes_prep:
    input:
        fai = REF_FASTA + ".fai",
        samples = expand(MAP_DIR + "{S}.sorted.nodup.nm.all.bam", S = ID)
    output:
        filter_fai = VCF_DIR_REF + REF_NAME + ".filter." + MIN_SIZE_SCAFFOLD + ".fasta.fai",
        regions = VCF_DIR_REF + REF_NAME + ".100kbp.regions"
    params:
        MIN_SIZE_SCAFFOLD
    threads: 1
    shell:
        """
	cat {input.fai} | awk '$2>= {params} {{print $0}}' > {output.filter_fai}
        code/split_ref_by_bai_datasize.py {input.samples} -r {input.fai} -s 100000 | sed 's/\t/:/' | sed 's/\t/-/' > {output.regions}
        """


rule freebayes_parallel:
    input:
        ref = REF_FASTA,
        regions = VCF_DIR_REF + REF_NAME + ".100kbp.regions",
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
        mkdir -p {params.tmpdir}
        export TMPDIR={params.tmpdir}
        freebayes-parallel {input.regions} {threads} -f {input.ref} {input.samples} --use-best-n-alleles 4 > {output.vcf}
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
        het = temp(VCF_DIR + SPECIES + ".heterozygosity.bed")
    log: VCF_DIR + SPECIES + ".proportion_heterozygosity.log"
    params:
        hetero = expand("het:{heterogametic}", heterogametic = HETEROGAMETIC),
        homo = expand("homo:{homogametic}", homogametic = HOMOGAMETIC)
    threads: 1
    shell:
        """
        python3 code/heterozygosity_per_indv.py {input} {output.het} {params.hetero} {params.homo} > {log}
        """

rule proportion_heterozygosity_window_old:
    input:
        het_sorted = VCF_DIR + SPECIES + ".heterozygosity.bed",
        windows = GENCOV_DIR_REF + "genome_5kb_windows.out"
    output:
        het_sorted_window = temp(VCF_DIR + SPECIES + ".heterozygosity.5kb.windows.bed"),
        het_sorted_window_mean = VCF_DIR + SPECIES + ".heterozygosity.5kb.windows.mean.bed",
	het_sexAverage = VCF_DIR + SPECIES + ".heterozygosity.sexAverage.bed"
    params:
        hetero = expand("het:{heterogametic}", heterogametic = HETEROGAMETIC),
        homo = expand("homo:{homogametic}", homogametic = HOMOGAMETIC)
    threads: 1
    shell:
        """
        bedtools intersect -a {input.windows} -b {input.het_sorted} -wa -wb | cut -f 1-3,7- > {output.het_sorted_window}
        cat {output.het_sorted_window} | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' | awk ' {{c[$1]++; for (i=2;i<=NF;i++) {{ s[$1"."i]+=$i}}; }} END {{for (k in c) {{printf "%s\t", k; for(i=2;i<NF;i++) printf "%.1f\\t", s[k"."i]/c[k]; printf "%.1f\\n", s[k"."NF]/c[k];}}}}' |  sed 's/STARTCOORD/\t/' | sed 's/ENDCOORD/\t/' | sed 's/ /\t/g' > {output.het_sorted_window_mean}
	python3 code/mean_heterozygosity_per_sex.py {output.het_sorted_window_mean} no-synteny {params.hetero} {params.homo} > {output.het_sexAverage}
	"""

rule proportion_heterozygosity_window:
    input:
        het_sorted = VCF_DIR + SPECIES + ".heterozygosity.bed",
        windows = GENCOV_DIR_REF + "genome_5kb_windows.out"
    output:
        het_sorted_window = temp(VCF_DIR + SPECIES + ".heterozygosity.5kb.windows.bed"),
        het_sorted_window_mean = VCF_DIR + SPECIES + ".heterozygosity.5kb.windows.NR.bed",
	het_sexAverage = VCF_DIR + SPECIES + ".heterozygosity.sexAverage.NR.bed"
    params:
        hetero = expand("het:{heterogametic}", heterogametic = HETEROGAMETIC),
        homo = expand("homo:{homogametic}", homogametic = HOMOGAMETIC)
    threads: 1
    shell:
        """
        bedtools intersect -a {input.windows} -b {input.het_sorted} -wa -wb | cut -f 1-3,7- > {output.het_sorted_window}
        cat {output.het_sorted_window} | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' | awk '
        BEGIN {{ FS=OFS="\t" }}
        {{for(i=2;i<=NF;i++)  
                a[$1][i]+=$i   
            a[$1][1]=i                    
        }}
        END {{
            for(i in a) {{
                for((j=2)&&b="";j<a[i][1];j++) 
                    b=b (b==""?"":OFS)a[i][j]
                print i,b          
            }}
        }} ' |  sed 's/STARTCOORD/\t/' | sed 's/ENDCOORD/\t/' | sed 's/ /\t/g' | awk '{{for(i=4;i<=NF;i++)$i/=50}}1' | sed 's/ /\t/g' > {output.het_sorted_window_mean}
	python3 code/mean_heterozygosity_per_sex.py {output.het_sorted_window_mean} no-synteny {params.hetero} {params.homo} > {output.het_sexAverage}
	"""

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
        out_scatter = RESULTDIR + SPECIES + "{synteny}{gencov}.scatter.{bp}bp.pdf",
        chromosomes = CHROMOSOMES,
	chromosomes_highlight = CHROMOSOMES_HIGHLIGHT,
	ED = expand("{ED}", ED = EDIT_DIST),
    wildcard_constraints:
        gencov = "gencov"
    shell:
        """
        Rscript code/plot_windows.R {input.cov} {input.snp} {params.out_scatter} {params.chromosomes} {params.chromosomes_highlight} {params.ED} 
        """

rule plotting_linear:
    input:
        cov = expand(RESULTDIR + SPECIES + "{{synteny}}{{gencov}}.nodup.nm.{ED}.{{bp}}bp.out", ED = EDIT_DIST),
        snp = RESULTDIR + SPECIES + "{synteny}diffHeterozygosity.{bp}bp.out"
    output:
        touch(RESULTDIR + SPECIES + "{synteny}{gencov}.plotting.linear.{bp}bp.done")
    threads: 1
    params:
        absolute_out = RESULTDIR + SPECIES + "{synteny}{gencov}.sexSpecificValues.{bp}bp.pdf",
        diff_out = RESULTDIR + SPECIES + "{synteny}{gencov}.sexDifference.{bp}bp.pdf",
        chromosomes = CHROMOSOMES,
	ED = expand("{ED}", ED = EDIT_DIST),
	nr_chromosomes = 50
    wildcard_constraints:
        gencov = "gencov"
    shell:
        """
        Rscript code/plot_windows_linear.R {input.cov} {input.snp} {params.absolute_out} {params.diff_out} {params.chromosomes} {params.ED} {params.nr_chromosomes}
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
        chromosomes = CHROMOSOMES,
	ED = expand("{ED}", ED = EDIT_DIST)
    threads: 1
    shell:
        """
        Rscript code/scatterplot_chr.R {input.cov} {input.snp} {params.out_scatter2D} {params.out_scatter3D} {params.chromosomes} {params.ED}
        """

##########################################################
################### MODIFY REF GENOME ####################
##########################################################

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
