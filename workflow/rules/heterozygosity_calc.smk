rule proportion_heterozygosity:
    input:
        outdir + "variant_calling/" + ref_genome_name_simple + ".biallelic.minQ20.minDP3.vcf.gz"
    output:
        het = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.bed")
    log: 
        outdir + "variant_calling/" + ref_genome_name_simple + ".proportion_heterozygosity.log"
    params:
        hetero = expand("het:{u.sample}__{u.unit}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.unit}", u=homogametic.itertuples())
    resources:
        mem_mb=2048
    shell:
        """
        python3 code/heterozygosity_per_indv.py {input} {output.het} {params.hetero} {params.homo} > {log}
        """


rule proportion_heterozygosity_window:
    input:
        het_sorted = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.bed",
        windows = outdir + "coverage/" + "genome_5kb_windows.out"
    output:
        het_sorted_window = temp(outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.bed"),
        het_sorted_window_mean = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.5kb.windows.NR.bed",
	    het_sexAverage = outdir + "variant_calling/" + ref_genome_name_simple + ".heterozygosity.sexAverage.NR.bed"
    params:
        hetero = expand("het:{u.sample}__{u.unit}", u=heterogametic.itertuples()),
        homo = expand("homo:{u.sample}__{u.unit}", u=homogametic.itertuples())
    resources:
        mem_mb=7048
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
