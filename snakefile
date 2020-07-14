from snakemake.utils import min_version
min_version("4.4.0")

# Written by Hanna Sigeman, Nov 2019

dir_path = os.getcwd()

###############################################
################## PATHS ######################
###############################################

ID = config["samples"]
SPECIES = config["species"]
HETEROGAMETIC = config["heterogametic"]
HOMOGAMETIC = config["homogametic"]
FQ_DIR = config["fastq"]
REF_SPECIES = config["ref_species"]
REF_DIR = config["ref_dir"]
REF_NAME = config["ref_name"]
PREFIX = SPECIES + "_ref_" + REF_SPECIES

REF_PATH = REF_DIR + REF_NAME
REF_FASTA = REF_DIR + REF_NAME + ".fasta"
MAP_DIR = "intermediate/bwa/" + REF_NAME + "/"
GENCOV_DIR = "intermediate/bedtools/" + PREFIX + "/"
GENCOV_DIR_REF = "intermediate/bedtools/" + REF_NAME + "/"
VCF_DIR = "intermediate/freebayes/" + PREFIX + "/"
VCF_DIR_REF = "intermediate/freebayes/" + REF_NAME + "/"
MATCHDIR = "intermediate/synteny_match/" + PREFIX + "/"
MATCHDIR_REF = "intermediate/synteny_match/" + REF_NAME + "/"
RESULTDIR = "results/" + PREFIX + "/"

SYNTENY_SPECIES = config["synteny_species"]
#COMP_GEN_SYNS = "intermediate/lastal_" + SYNTENY_SPECIES + "/" + PREFIX + "/"
COMP_GEN_SYNS = "intermediate/lastal_" + SYNTENY_SPECIES + "/" + REF_NAME + "/"
SYNS_DB = "data/meta/my" + SYNTENY_SPECIES + "db"

EDIT_DIST = ["all", "0.0", "0.1", "0.2", "0.3", "0.4"]

###############################################
################## RULES ######################
###############################################

rule all: 
     input: 
        REF_FASTA + ".bwt",
        REF_FASTA + ".fai" ,
        expand(MAP_DIR + "{S}.sorted.status", S = ID),
        expand(MAP_DIR + "{S}.sorted.nodup.nm.all.status", S = ID),
        expand(MAP_DIR + "{S}.sorted.nodup.nm.{ED}.bam.bai", S = ID, ED = EDIT_DIST),
        VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf",
        COMP_GEN_SYNS + REF_SPECIES + "_align_converted",
        expand(MAP_DIR + "{S}.sorted.flagstat", S = ID),
        expand(MAP_DIR + "{S}.sorted.nodup.nm.{ED}.flagstat", S = ID, ED = EDIT_DIST),
        MATCHDIR_REF + "genome_windows.out",
        MATCHDIR_REF + "bestMatch.status",
        expand(MATCHDIR + "gencov.nodup.nm.{ED}.norm." + SYNTENY_SPECIES + ".out", ED = EDIT_DIST),
        #REF_DIR + REF_NAME + "_nonRefAc_consensus.fasta",
        RESULTDIR + SPECIES + ".circlize.pdf",
        VCF_DIR + SPECIES + ".diffHeterozygosity.bed"

##########################################################  
##################### INDEX GENOME #######################      
##########################################################

rule index_fasta_bwa:
    input: 
        ref = REF_FASTA
    output:
        ref_bwt = REF_FASTA + ".bwt" 
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
        ref = REF_FASTA
    output: 
        ref_fai = REF_FASTA + ".fai"
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
        out = temp(MAP_DIR + "{S}.sorted.bam"),
        log = MAP_DIR + "{S}.sorted.status"
    threads: 3
    params:
        tmpdir = MAP_DIR + "{S}_temp_sort/"
    shell:
        """
        mkdir {params.tmpdir}
       samtools sort -@ {threads} {input} -T {params.tmpdir} > {output.out}
        rm -r {params.tmpdir}
        echo "DONE" > {output.log}
        """

rule remove_duplicates: 
    input: 
        MAP_DIR + "{S}.sorted.bam"
    output: 
        out = MAP_DIR + "{S}.sorted.nodup.nm.all.bam",
        log = MAP_DIR + "{S}.sorted.nodup.nm.all.status"
    params:
        tmpdir = MAP_DIR + "{S}_temp_dupl/"
    shell: 
        """
        mkdir {params.tmpdir}
        picard MarkDuplicates -Xmx10g MAX_FILE_HANDLES=500 REMOVE_DUPLICATES=true I={input} O={output.out} M={input}_duplicatedata.txt TMP_DIR={params.tmpdir}
        rm -r {params.tmpdir}
        echo "DONE" > {output.log}
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
      MAP_DIR + "{S}.sorted{V}bam" 		# V = ".", ".nodup.nm.all.", ".nodup.nm.0.0.", ".nodup.nm.0.1.", ".nodup.nm.0.2.", ".nodup.nm.0.3.", ".nodup.nm.0.4."
    output:
      MAP_DIR + "{S}.sorted{V}flagstat",
    shell:
        """
        samtools flagstat {input} > {output}  
        """

##########################################################  
#################### GENOME COVERAGE #####################       
########################################################## 

rule gencov_prepare_fasta:
    input: 
        ref_fai = REF_FASTA + ".fai"
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
        GENCOV_DIR + "gencov.nodup.nm.{V}.out"
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
        f = expand(MAP_DIR + "{heterogametic}.sorted.nodup.nm.all.bam", heterogametic = HETEROGAMETIC),
        m = expand(MAP_DIR + "{homogametic}.sorted.nodup.nm.all.bam", homogametic = HOMOGAMETIC)
    output: 
        vcf = VCF_DIR + SPECIES + ".vcf",
        log = VCF_DIR + SPECIES + ".vcf.status"
    threads: 18
    params:
        tmpdir=VCF_DIR + "temp/"
    shell: 
        """
        mkdir {params.tmpdir}
        export TMPDIR={params.tmpdir}
        freebayes-parallel {input.regions} {threads} -f {input.ref} {input.f} {input.m} > {output.vcf}
        rm -r {params.tmpdir}
        echo "DONE" > {output.log}
        """

rule vcftools_filter:
     input:
         VCF_DIR + SPECIES + ".vcf"
     output:
         VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf"
     threads: 1
     shell:
         """
         vcftools --vcf {input} --min-alleles 2 --max-alleles 2 --remove-filtered-geno-all --minQ 20 --minDP 3 --recode --stdout > {output}
         """

rule proportion_heterozygosity:
     input:
         vcf = VCF_DIR + SPECIES + ".biallelic.minQ20.minDP3.vcf"
     output:
         VCF_DIR + SPECIES + ".diffHeterozygosity.bed"
     params:
         hetero = expand("{heterogametic}", heterogametic = HETEROGAMETIC),
         homo = expand("{homogametic}", homogametic = HOMOGAMETIC)
     threads: 1
     shell:
         """
         python3 code/calculate_prop_heterozygosity.py {input.vcf} {output} {params.hetero} {params.homo}
         """

##########################################################  
#################### SYNTENY ANALYSIS ####################      
########################################################## 

rule fasta_formatter:
    input: 
        ref = REF_FASTA
    output: 
        REF_PATH + "wrap.fasta"
    threads: 1
    shell: 
        """
        fasta_formatter -i {input} -w 80 -o {output}
        """

rule lastal_syns:
    input: 
        REF_PATH + "wrap.fasta"
    output: 
        COMP_GEN_SYNS + REF_SPECIES + "_align"
    params: 
        db = SYNS_DB
    threads: 15
    shell: 
        """
        lastal {params.db} {input} -P 15 | last-split > {output}
        """

rule maf_convert_syns:
    input: 
        COMP_GEN_SYNS + REF_SPECIES + "_align"
    output: 
        COMP_GEN_SYNS + REF_SPECIES + "_align_converted"
    threads: 1
    shell: 
        """
        maf-convert psl {input} > {output}
        """

##########################################################  
################### MATCHING DATASETS ####################      
########################################################## 

rule matchScaffold2Chr:
    input:
        syns = COMP_GEN_SYNS + REF_SPECIES + "_align_converted",
        gencov = GENCOV_DIR_REF + "genome_5kb_windows.out"
    output:
        windows = MATCHDIR_REF  + "genome_windows.out",
        bestMatch = MATCHDIR_REF + "bestMatch.list",
        log = MATCHDIR_REF + "bestMatch.status"
    params:
        temp = MATCHDIR_REF + "temp/",
        abswindow = dir_path + "/" + MATCHDIR_REF  + "genome_windows.out",
        absBestMatch = dir_path + "/" + MATCHDIR_REF + "bestMatch.list",
        windowsfile= "genome_windows.out",
        absLog = dir_path + "/" + MATCHDIR_REF + "bestMatch.status"
    shell:
        """
        cat {input.syns} | awk '{{print $10,$12,$13,$14,$16,$17,$1}}' | sed 's/ /\t/g' | bedtools intersect -a stdin -b {input.gencov} -wa -wb | awk '{{if($10-$9==\"5000\") print $8,$9,$10,$7,$1,$2,$3,$4,$5,$6}}' | sed 's/ /\t/g' | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' > {output.windows}

        mkdir {params.temp}
        cd {params.temp}
        
        awk \'{{print >> $1; close($1)}}\' {params.abswindow}

        ls | while read file; do cat $file | awk -v max=0 '{{if($2>max){{want=$0; max=$2}}}}END{{print want}}' ; done | sed 's/STARTCOORD/\t/' | sed 's/ENDCOORD/\t/' >> {params.absBestMatch} 

        cd {dir_path}

        rm -r {params.temp}

        echo "DONE" > {params.absLog}
        """

rule matchScaffold2Chr_snp:
    input:
        bestMatch = MATCHDIR_REF + "bestMatch.list",
        heterozygosity = VCF_DIR + SPECIES + ".diffHeterozygosity.bed"
    output:
        bestmatch = MATCHDIR + SPECIES + ".diffHeterozygosity.bestMatch." + SYNTENY_SPECIES,
        bestmatch_small = MATCHDIR + SPECIES + ".diffHeterozygosity.bestMatch." + SYNTENY_SPECIES + ".small"
    shell:
        """
        bedtools intersect -a {input.heterozygosity} -b {input.bestMatch} -wa -wb > {output.bestmatch}
        
        cat {output.bestmatch} | cut -f 4,12,13,14 > {output.bestmatch_small}
        """

rule matchScaffold2Chr_cov:
    input:
        bestMatch = MATCHDIR_REF + "bestMatch.list",
        cov = GENCOV_DIR + "gencov.nodup.nm.{ED}.out" # ED = all, 0.0, 0.1, 0.2, 0.3, 0.4
    output:
        bestMatch = MATCHDIR + "gencov.nodup.nm.{ED}." + SYNTENY_SPECIES + ".out",
        bestMatch_norm = MATCHDIR + "gencov.nodup.nm.{ED}.norm." + SYNTENY_SPECIES + ".out"
    threads: 1
    shell:
        """
        bedtools intersect -a {input.bestMatch} -b {input.cov} -wa -wb > {output.bestMatch}

        python3 code/matchScaffold2chr_cov.py {output.bestMatch} > {output.bestMatch_norm}
        """

rule confirm_sexing:
    input:
        MATCHDIR + "gencov.nodup.nm.all." + SYNTENY_SPECIES + ".out"
    output:
        RESULTDIR + SPECIES + ".gencovIndv.pdf"
    threads: 1
    shell:
        """
        Rscript ...
        """

##########################################################  
################### MODIFY REF GENOME ####################      
########################################################## 

# Use major allele instead or HISAT2
#rule modify_genome:
#    input: 
#        vcf = VCF_DIR + SPECIES + ".vcf",
#        ref = REF_FASTA
#    output: 
#        vcf = VCF_DIR + SPECIES + ".non-ref-ac_2_biallelic_qual.vcf",
#        gz = VCF_DIR + SPECIES + ".non-ref-ac_2_biallelic_qual.vcf.gz",
#        ref = REF_DIR + REF_NAME + "_nonRefAc_consensus.fasta"
#    threads: 1
#    shell: 
#        """
#        vcftools --vcf {input.vcf} --non-ref-ac 2 --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout > {output.vcf}
#        bgzip -c {output.vcf} > {output.gz}
#        tabix -p vcf {output.gz}
#        cat {input.ref} | bcftools consensus {output.gz} > {output.ref}
#        """


###########################################################
################# STATISTICAL CALCULATIONS ##################
###########################################################

rule calculate_heterozygosity:
    input:
        MATCHDIR + SPECIES + ".diffHeterozygosity.bestMatch." + SYNTENY_SPECIES + ".small"
    output:
        Mb = RESULTDIR + SPECIES + ".diffHeterozygosity.1Mbp.out",
        kb = RESULTDIR + SPECIES + ".diffHeterozygosity.100kbp.out"
    threads: 1
    shell:
        """
        Rscript code/calculate_heterozygosityDiff_windows.R {input} {output.Mb} {output.kb}
        """

rule calculate_ratio:
    input:
        MATCHDIR + "{type}." + SYNTENY_SPECIES + ".out"
    output:
        Mb = RESULTDIR + SPECIES + ".{type}.1Mbp.out",
        kb = RESULTDIR + SPECIES + ".{type}.100kbp.out"
    threads: 1
    shell:
        """
        Rscript code/calculate_gencov_windows.R {input} {output.Mb} {output.kb}
        """

rule plotting:
    input: 
        cov0 = RESULTDIR + SPECIES + ".gencov.nodup.nm.0.0.norm.1Mbp.out",
        cov2 = RESULTDIR + SPECIES + ".gencov.nodup.nm.0.2.norm.1Mbp.out",
        cov4 = RESULTDIR + SPECIES + ".gencov.nodup.nm.0.4.norm.1Mbp.out",
        snp = RESULTDIR + SPECIES + ".diffHeterozygosity.1Mbp.out"
    output: 
        RESULTDIR + SPECIES + ".circlize.pdf"
    threads: 1
    shell: 
        """
        Rscript code/circlize_plot.R {input.cov0} {input.cov2} {input.cov4} {input.snp} {output}
        """


