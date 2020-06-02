from snakemake.utils import min_version
min_version("4.4.0")

# Written by Hanna Sigeman, Nov 2019

dir_path = os.getcwd()

###############################################
################## PATHS ######################
###############################################

ID = config["samples"]
SPECIES = config["species"]
FEMALE = config["female"]
MALE = config["male"]
FQ_DIR = config["fastq"]
REF_SPECIES = config["ref_species"]
REF_DIR = config["ref_dir"]
REF_NAME = config["ref_name"]
PREFIX = SPECIES + "_ref_" + REF_SPECIES
NM_VALUES = [0, 1, 2, 3, 4]

REF_PATH = REF_DIR + REF_NAME
REF_FASTA = REF_DIR + REF_NAME + ".fasta"
MAP_DIR = "intermediate/bwa/" + PREFIX + "/" 
GENCOV_DIR = "intermediate/bedtools/" + PREFIX + "/"
VCF_DIR = "intermediate/freebayes/" + PREFIX + "/"
MATCHDIR = "intermediate/synteny_match/" + PREFIX + "/"
RESULTDIR = "results/" + PREFIX + "/"

COMP_GEN_ZF = "intermediate/lastal_ZF/" + PREFIX + "/"
ZF_DB = "data/meta/myZFdb"

###############################################
################## RULES ######################
###############################################

rule all: 
     input: 
        REF_FASTA + ".bwt",
        REF_FASTA + ".fai" ,
        expand(MAP_DIR + "{S}" + ".sorted.status", S = ID),
        expand(MAP_DIR + "{S}" + ".sorted.nodup.status", S = ID),
        expand(MAP_DIR + "{S}" + ".sorted.nodup.bam.bai", S = ID),
        expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0." + "{NM}" + ".bam.bai", S = ID, NM = NM_VALUES),
        #expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0.1.bam.bai", S = ID),
        #expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0.2.bam.bai", S = ID),
        #expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0.3.bam.bai", S = ID),
        #expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0.4.bam.bai", S = ID),
        VCF_DIR + SPECIES + ".vcf.status",
        COMP_GEN_ZF + SPECIES + "_align_converted",
        expand(MAP_DIR + "{S}" + ".sorted.flagstat", S = ID),
        expand(MAP_DIR + "{S}" + ".sorted.nodup.flagstat", S = ID),
        expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0." + "{NM}" + ".flagstat", S = ID, NM = NM_VALUES),
        #expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0.1.flagstat", S = ID),
        #expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0.2.flagstat", S = ID),
        #expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0.3.flagstat", S = ID),
        #expand(MAP_DIR + "{S}" + ".sorted.nodup.nm.0.4.flagstat", S = ID),
        MATCHDIR + "genome_windows.out",
        MATCHDIR + "bestMatch.status",
        MATCHDIR + SPECIES + ".singleton.bestMatch.zf.small",
        MATCHDIR + "gencov.nodup.nm.all.zf.out",
        expand(MATCHDIR + "gencov.nodup.nm.0." + "{NM}" + ".zf.out", NM = NM_VALUES),
        #MATCHDIR + "gencov.nodup.nm.0.1.zf.out",
        #MATCHDIR + "gencov.nodup.nm.0.2.zf.out",
        #MATCHDIR + "gencov.nodup.nm.0.3.zf.out",
        #MATCHDIR + "gencov.nodup.nm.0.4.zf.out",
        VCF_DIR + SPECIES + ".non-ref-ac_2_biallelic_qual.vcf",
        VCF_DIR + SPECIES + ".non-ref-ac_2_biallelic_qual.vcf.gz",
        REF_DIR + REF_NAME + "_nonRefAc_consensus.fasta",
        RESULTDIR + SPECIES + ".circlize.pdf"


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

# Is not used
#rule count_reads: 
#    input: 
#        R1= FQ_DIR + "{S}_reverse_paired.fq.gz"
#    output: 
#        GENCOV_DIR + "{S}" + "_fastq_readcount.out"
#    shell:
#        """ 
#        echo {input.R1} | tr '\\n' '\\t' ; zcat {input.R1} | echo $((`wc -l`/4)) > {output}
#        """ 


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
        temp(MAP_DIR + "{S}" + ".bam")
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
        MAP_DIR + "{S}" + ".bam"
    output:
        out = temp(MAP_DIR + "{S}" + ".sorted.bam"),
        log = MAP_DIR + "{S}" + ".sorted.status"
    threads: 3
    params:
        tmpdir = MAP_DIR + "{S}" + "_temp_sort/"
    shell:
        """
        mkdir {params.tmpdir}
       samtools sort -@ {threads} {input} -T {params.tmpdir} > {output.out}
        rm -r {params.tmpdir}
        echo "DONE" > {output.log}
        """

rule remove_duplicates: 
    input: 
        MAP_DIR + "{S}" + ".sorted.bam"
    output: 
        out = MAP_DIR + "{S}" + ".sorted.nodup.bam",
        log = MAP_DIR + "{S}" + ".sorted.nodup.status"
    params:
        tmpdir = MAP_DIR + "{S}" + "_temp_dupl/"
    shell: 
        """
        mkdir {params.tmpdir}
        picard MarkDuplicates -Xmx10g MAX_FILE_HANDLES=500 REMOVE_DUPLICATES=true I={input} O={output.out} M={input}_duplicatedata.txt TMP_DIR={params.tmpdir}
        rm -r {params.tmpdir}
        echo "DONE" > {output.log}
        """

rule index_bam: 
    input: 
        MAP_DIR + "{S}" + ".sorted.nodup.bam"
    output: 
        MAP_DIR + "{S}" + ".sorted.nodup.bam.bai"
    threads: 1
    shell: 
        """
        samtools index {input}
        """

# not used
#rule count_reads_samtools: 
#    input: 
#        MAP_DIR + "{S}" + ".sorted.nodup.bam"
#    output: 
#        GENCOV_DIR + "{S}" + "_fastq_readcount_samtools.out"
#    shell:
#        """ 
#        samtools view -c -F 4 {input} > {output}
#        """ 

rule mismatch_bam: 
    input: 
        bam = MAP_DIR + "{S}" + ".sorted.nodup.bam", 
        ref = REF_FASTA
    output: 
        MAP_DIR + "{S}" + ".sorted.nodup.nm.0." + "{NM}" + ".bam"
    threads: 2
    params:
        "\"NM:i:[0-{NM}]\"" #the NM values need to get inserted in this params
    shell: 
        """
        samtools view -@ {threads} {input.bam} | grep -E {params} | samtools view -bS -T {input.ref} - > {output}
        """

rule index_mismatch_bam: 
    input: 
        MAP_DIR + "{S}" + ".sorted.nodup.nm.0." + "{NM}" + ".bam"
    output: 
        MAP_DIR + "{S}" + ".sorted.nodup.nm.0." + "{NM}" + ".bam.bai"
    threads: 1
    shell: 
        """
        samtools index {input}
        """


#rule mismatch_bam1: 
#    input: 
#        bam = MAP_DIR + "{S}" + ".sorted.nodup.bam", 
#        ref = REF_FASTA
#    output: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.1.bam"
#    threads: 2
#    params:
#        "\"NM:i:[0-1]\""
#    shell: 
#        """
#        samtools view -@ {threads} {input.bam} | grep -E {params} | samtools view -bS -T {input.ref} - > {output}
#        """

#rule index_mismatch_bam1: 
#    input: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.1.bam"
#    output: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.1.bam.bai"
#    threads: 1
#    shell: 
#        """
#        samtools index {input}
#        """

#rule mismatch_bam2: 
#    input: 
#        bam = MAP_DIR + "{S}" + ".sorted.nodup.bam", 
#        ref = REF_FASTA
#    output: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.2.bam"
#    threads: 2
#    params:
#        "\"NM:i:[0-2]\""
#    shell: 
#        """
#        samtools view -@ {threads} {input.bam} | grep -E {params} | samtools view -bS -T {input.ref} - > {output}
#        """

#rule index_mismatch_bam2: 
#    input: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.2.bam"
#    output: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.2.bam.bai"
#    threads: 1
#    shell: 
#        """
#        samtools index {input}
#        """

#rule mismatch_bam3: 
#    input: 
#        bam = MAP_DIR + "{S}" + ".sorted.nodup.bam", 
#        ref = REF_FASTA
#    output: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.3.bam"
#    threads: 2
#    params:
#        "\"NM:i:[0-3]\""
#    shell: 
#        """
#        samtools view -@ {threads} {input.bam} | grep -E {params} | samtools view -bS -T {input.ref} - > {output}
#        """

#rule index_mismatch_bam3: 
#    input: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.3.bam"
#    output: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.3.bam.bai"
#    threads: 1
#    shell: 
#        """
#        samtools index {input}
#        """

#rule mismatch_bam4: 
#    input: 
#        bam = MAP_DIR + "{S}" + ".sorted.nodup.bam", 
#        ref = REF_FASTA
#    output: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.4.bam"
#    threads: 2
#    params:
#        "\"NM:i:[0-4]\""
#    shell: 
#        """
#        samtools view -@ {threads} {input.bam} | grep -E {params} | samtools view -bS -T {input.ref} - > {output}
#        """

#rule index_mismatch_bam4: 
#    input: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.4.bam"
#    output: 
#        MAP_DIR + "{S}" + ".sorted.nodup.nm.0.4.bam.bai"
#    threads: 1
#    shell: 
#        """
#        samtools index {input}
#        """

rule flagstat:
    input:
      expand(MAP_DIR + "{S}" + ".sorted" + "{VAL}" + ".bam", S = ID, VAL = ["", ".nodup", ".nodup.nm.0.0", ".nodup.nm.0.1",".nodup.nm.0.2",".nodup.nm.0.3",".nodup.nm.0.4"])
      #sort = MAP_DIR + "{S}" + ".sorted.bam",
      #nodup = MAP_DIR + "{S}" + ".sorted.nodup.bam",
      #zero =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.0.bam",
      #one =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.1.bam",
      #two =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.2.bam",
      #three =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.3.bam",
      #four =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.4.bam"
    output:
      expand(MAP_DIR + "{S}" + ".sorted" + "{VAL}" + ".flagstat", S = ID, VAL = ["", ".nodup", ".nodup.nm.0.0", ".nodup.nm.0.1",".nodup.nm.0.2",".nodup.nm.0.3",".nodup.nm.0.4"])
      #sort = MAP_DIR + "{S}" + ".sorted.flagstat",
      #nodup = MAP_DIR + "{S}" + ".sorted.nodup.flagstat",
      #zero =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.0.flagstat",
      #one =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.1.flagstat",
      #two =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.2.flagstat",
      #three =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.3.flagstat",
      #four =  MAP_DIR + "{S}" + ".sorted.nodup.nm.0.4.flagstat"
    shell:
        """
        samtools flagstat {input} > {output} 
        #samtools flagstat {input.sort} > {output.sort}
        #samtools flagstat {input.nodup} > {output.nodup} 
        #samtools flagstat {input.zero} > {output.zero} 
        #samtools flagstat {input.one} > {output.one} 
        #samtools flagstat {input.two} > {output.two} 
        #samtools flagstat {input.three} > {output.three} 
        #samtools flagstat {input.four} > {output.four} 
        """

##########################################################  
#################### GENOME COVERAGE #####################       
########################################################## 

rule gencov_prepare_fasta:
    input: 
        ref_fai = REF_FASTA + ".fai"
    output: 
        GENCOV_DIR + "genome_5kb_windows.out"
    threads: 1
    shell: 
        """
        bedtools makewindows -g {input} -w 5000 -s 5000 > {output}
        """

rule gencov_bedtoolsall:
    input: 
        bam = expand(MAP_DIR + "{S}" + ".sorted.nodup.bam", S = ID),
        bai = expand(MAP_DIR + "{S}" + ".sorted.nodup.bam.bai", S = ID),
#        bam_f = MAP_DIR + FEMALE + ".sorted.nodup.bam",
#        bai_f = MAP_DIR + FEMALE + ".sorted.nodup.bam.bai",
#        bam_m = MAP_DIR + MALE + ".sorted.nodup.bam",
#        bai_m = MAP_DIR + MALE + ".sorted.nodup.bam.bai",
        bed = GENCOV_DIR + "genome_5kb_windows.out"
    output: 
        GENCOV_DIR + "gencov.nodup.nm.all.out"
    threads: 2
    shell: 
        """
        bedtools multicov -bams {input.bam} -bed {input.bed} > {output}
#        bedtools multicov -bams {input.bam_f} {input.bam_m} -bed {input.bed} > {output}
        """

rule gencov_bedtools0:
    input: 
        bam_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.0.bam",
        bai_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.0.bam.bai",
        bam_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.0.bam",
        bai_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.0.bam.bai",
        bed = GENCOV_DIR + "genome_5kb_windows.out"
    output: 
        GENCOV_DIR + "gencov.nodup.nm.0.0.out"
    threads: 2
    shell: 
        """
        bedtools multicov -bams {input.bam_f} {input.bam_m} -bed {input.bed} > {output}
        """
   
rule gencov_bedtools1:
    input: 
        bam_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.1.bam",
        bai_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.1.bam.bai",
        bam_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.1.bam",
        bai_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.1.bam.bai",
        bed = GENCOV_DIR + "genome_5kb_windows.out"
    output: 
        GENCOV_DIR + "gencov.nodup.nm.0.1.out"
    threads: 2
    shell: 
        """
        bedtools multicov -bams {input.bam_f} {input.bam_m} -bed {input.bed} > {output}
        """

rule gencov_bedtools2:
    input: 
        bam_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.2.bam",
        bai_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.2.bam.bai",
        bam_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.2.bam",
        bai_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.2.bam.bai",
        bed = GENCOV_DIR + "genome_5kb_windows.out"
    output: 
        GENCOV_DIR + "gencov.nodup.nm.0.2.out"
    threads: 2
    shell: 
        """
        bedtools multicov -bams {input.bam_f} {input.bam_m} -bed {input.bed} > {output}
        """

rule gencov_bedtools3:
    input: 
        bam_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.3.bam",
        bai_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.3.bam.bai",
        bam_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.3.bam",
        bai_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.3.bam.bai",
        bed = GENCOV_DIR + "genome_5kb_windows.out"
    output: 
        GENCOV_DIR + "gencov.nodup.nm.0.3.out"
    threads: 2
    shell: 
        """
        bedtools multicov -bams {input.bam_f} {input.bam_m} -bed {input.bed} > {output}
        """

rule gencov_bedtools4:
    input: 
        bam_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.4.bam",
        bai_f = MAP_DIR + FEMALE + ".sorted.nodup.nm.0.4.bam.bai",
        bam_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.4.bam",
        bai_m = MAP_DIR + MALE + ".sorted.nodup.nm.0.4.bam.bai",
        bed = GENCOV_DIR + "genome_5kb_windows.out"
    output: 
        GENCOV_DIR + "gencov.nodup.nm.0.4.out"
    threads: 2
    shell: 
        """
        bedtools multicov -bams {input.bam_f} {input.bam_m} -bed {input.bed} > {output}
        """


##########################################################  
#################### VARIANT CALLING #####################      
########################################################## 


rule freebayes_prep:
    input: 
        ref_fai = REF_FASTA + ".fai"
    output: 
        VCF_DIR + SPECIES + ".100kbp.regions"
    threads: 4
    shell: 
        """
        fasta_generate_regions.py {input} 100000 > {output}
        """

rule freebayes_parallel:
    input: 
        ref = REF_FASTA,
        regions = VCF_DIR + SPECIES + ".100kbp.regions",
        f = MAP_DIR + FEMALE + ".sorted.nodup.bam",
        m = MAP_DIR + MALE + ".sorted.nodup.bam"
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

rule vcftools_singletons:
    input: 
        VCF_DIR + SPECIES + ".vcf"
    output: 
        VCF_DIR + SPECIES + ".singletons.bed"
    threads: 1
    shell: 
        """
        vcftools --vcf {input} --singletons --remove-filtered-geno-all --minQ 20 --minDP 3 --stdout | awk -v OFS=\"\\t\" \'{{print $1,$2,$2+1,$3,$4,$5}}\' > {output}
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

rule lastal_zf:
    input: 
        REF_PATH + "wrap.fasta"
    output: 
        COMP_GEN_ZF + SPECIES + "_align"
    params: 
        db = ZF_DB
    threads: 15
    shell: 
        """
        lastal {params.db} {input} -P 15 | last-split > {output}
        """

rule maf_convert_zf:
    input: 
        COMP_GEN_ZF + SPECIES + "_align"
    output: 
        COMP_GEN_ZF + SPECIES + "_align_converted"
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
        zf = COMP_GEN_ZF + SPECIES + "_align_converted",
        gencov = GENCOV_DIR + "genome_5kb_windows.out"
    output:
        windows = MATCHDIR  + "genome_windows.out",
        bestMatch = MATCHDIR + "bestMatch.list",
        log = MATCHDIR + "bestMatch.status"
    params:
        temp = MATCHDIR + "temp/",
        abswindow = dir_path + "/" + MATCHDIR  + "genome_windows.out",
        absBestMatch = dir_path + "/" + MATCHDIR + "bestMatch.list",
        windowsfile= "genome_windows.out",
        absLog = dir_path + "/" + MATCHDIR + "bestMatch.status"
    shell:
        """
        cat {input.zf} | awk '{{print $10,$12,$13,$14,$16,$17,$1}}' | sed 's/ /\t/g' | bedtools intersect -a stdin -b {input.gencov} -wa -wb | awk '{{if($10-$9==\"5000\") print $8,$9,$10,$7,$1,$2,$3,$4,$5,$6}}' | sed 's/ /\t/g' | sed 's/\t/STARTCOORD/' | sed 's/\t/ENDCOORD/' > {output.windows}

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
        bestMatch = MATCHDIR + "bestMatch.list",
        vcf = VCF_DIR + SPECIES + ".singletons.bed"
    output:
        bestMatch_singleton = MATCHDIR + SPECIES + ".singleton.bestMatch.zf",
        bestMatch_singleton_small = MATCHDIR + SPECIES + ".singleton.bestMatch.zf.small"
    shell:
        """
        cat {input.vcf} | grep -v CHROM | sed 's/ /\t/g' | bedtools intersect -a stdin -b {input.bestMatch} -wa -wb > {output.bestMatch_singleton}

        cat {output.bestMatch_singleton} | cut -f 4,5,6,14,15,16 > {output.bestMatch_singleton_small}
        """

rule matchScaffold2Chr_cov:
    input:
        bestMatch = MATCHDIR + "bestMatch.list",
        covall = GENCOV_DIR + "gencov.nodup.nm.all.out",
        cov0 = GENCOV_DIR + "gencov.nodup.nm.0.0.out",
        cov1 = GENCOV_DIR + "gencov.nodup.nm.0.1.out",
        cov2 = GENCOV_DIR + "gencov.nodup.nm.0.2.out",
        cov3 = GENCOV_DIR + "gencov.nodup.nm.0.3.out",
        cov4 = GENCOV_DIR + "gencov.nodup.nm.0.4.out"
    output:
        outcovall = MATCHDIR + "gencov.nodup.nm.all.zf.out",
        outcov0 = MATCHDIR + "gencov.nodup.nm.0.0.zf.out",
        outcov1 = MATCHDIR + "gencov.nodup.nm.0.1.zf.out",
        outcov2 = MATCHDIR + "gencov.nodup.nm.0.2.zf.out",
        outcov3 = MATCHDIR + "gencov.nodup.nm.0.3.zf.out",
        outcov4 = MATCHDIR + "gencov.nodup.nm.0.4.zf.out"
    shell:
        """
        join -j 1 -o 1.1,1.2,1.7,1.8,1.9,2.3,2.4 <(cat {input.bestMatch} | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) <(cat {input.covall} | awk '$3-$2=="5000" {{ print $0}}' | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) | sed 's/STARTCOORD/\t/' | sed 's/\ /\t/g' > {output.outcovall}
        join -j 1 -o 1.1,1.2,1.7,1.8,1.9,2.3,2.4 <(cat {input.bestMatch} | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) <(cat {input.cov0} | awk '$3-$2=="5000" {{ print $0}}' | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) | sed 's/STARTCOORD/\t/' | sed 's/\ /\t/g' > {output.outcov0}
        join -j 1 -o 1.1,1.2,1.7,1.8,1.9,2.3,2.4 <(cat {input.bestMatch} | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) <(cat {input.cov1} | awk '$3-$2=="5000" {{ print $0}}' | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) | sed 's/STARTCOORD/\t/' | sed 's/\ /\t/g' > {output.outcov1}
        join -j 1 -o 1.1,1.2,1.7,1.8,1.9,2.3,2.4 <(cat {input.bestMatch} | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) <(cat {input.cov2} | awk '$3-$2=="5000" {{ print $0}}' | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) | sed 's/STARTCOORD/\t/' | sed 's/\ /\t/g' > {output.outcov2}
        join -j 1 -o 1.1,1.2,1.7,1.8,1.9,2.3,2.4 <(cat {input.bestMatch} | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) <(cat {input.cov3} | awk '$3-$2=="5000" {{ print $0}}' | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) | sed 's/STARTCOORD/\t/' | sed 's/\ /\t/g' > {output.outcov3}
        join -j 1 -o 1.1,1.2,1.7,1.8,1.9,2.3,2.4 <(cat {input.bestMatch} | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) <(cat {input.cov4} | awk '$3-$2=="5000" {{ print $0}}' | sed 's/\t/STARTCOORD/' | sort -k 1b,1 -k2,2) | sed 's/STARTCOORD/\t/' | sed 's/\ /\t/g' > {output.outcov4}
        """

##########################################################  
################### MODIFY REF GENOME ####################      
########################################################## 

rule modify_genome:
    input: 
        vcf = VCF_DIR + SPECIES + ".vcf",
        ref = REF_FASTA
    output: 
        vcf = VCF_DIR + SPECIES + ".non-ref-ac_2_biallelic_qual.vcf",
        gz = VCF_DIR + SPECIES + ".non-ref-ac_2_biallelic_qual.vcf.gz",
        ref = REF_DIR + REF_NAME + "_nonRefAc_consensus.fasta"
    threads: 1
    shell: 
        """
        vcftools --vcf {input.vcf} --non-ref-ac 2 --min-alleles 2 --max-alleles 2 --remove-filtered-all --recode --stdout > {output.vcf}
        bgzip -c {output.vcf} > {output.gz}
        tabix -p vcf {output.gz}
        cat {input.ref} | bcftools consensus {output.gz} > {output.ref}
        """

rule plotting:
    input: 
        cov0 = MATCHDIR + "gencov.nodup.nm.0.0.zf.out",
        cov2 = MATCHDIR + "gencov.nodup.nm.0.2.zf.out",
        cov4 = MATCHDIR + "gencov.nodup.nm.0.4.zf.out",
        snp = MATCHDIR + SPECIES + ".singleton.bestMatch.zf.small"
    output: 
        RESULTDIR + SPECIES + ".circlize.pdf"
    params:
        female = FEMALE,
        male = MALE,
        species = PREFIX
    threads: 1
    shell: 
        """
        Rscript code/Rscript_norm_plot.R {params.species} {params.female} {params.male} {input.cov0} {input.cov2} {input.cov4} {input.snp} {output}
        """
