#!/bin/bash -l
#SBATCH -A snic2019-3-374
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J all_freq

# Generates a file from which the allele frequency can be plotted.


module load bioinfo-tools vcftools/0.1.16 BEDTools/2.29.2 python3/3.7.2

# Allele frequency for each sex
#vcftools --vcf intermediate/freebayes/AlaArvIT_ref_AlaArv/AlaArvIT.biallelic.minQ20.minDP3.vcf --indv SJ-2333-Aarv-249-IT_S21_L002 --indv SJ-2333-Aarv-441-IT_S25_L002 --indv SJ-2333-Aarv-460-IT_S27_L002 --freq --stdout > intermediate/freebayes/AlaArvIT_ref_AlaArv/allFreq.heterogametic.out
#vcftools --vcf intermediate/freebayes/AlaArvIT_ref_AlaArv/AlaArvIT.biallelic.minQ20.minDP3.vcf --indv SJ-2333-Aarv-389-IT_S17_L002 --indv SJ-2333-Aarv-449-IT_S19_L002 --indv SJ-2333-Aarv-939-IT_S20_L002 --freq --stdout > intermediate/freebayes/AlaArvIT_ref_AlaArv/allFreq.homogametic.out

# Keep only sites with a certain allele frequency
python3 code/filter_allFreq.py intermediate/freebayes/AlaArvIT_ref_AlaArv/allFreq.homogametic.out  homogametic > intermediate/freebayes/AlaArvIT_ref_AlaArv/allFreq.out
python3 code/filter_allFreq.py intermediate/freebayes/AlaArvIT_ref_AlaArv/allFreq.heterogametic.out heterogametic >> intermediate/freebayes/AlaArvIT_ref_AlaArv/allFreq.out

# Matches skylark scaffolds too zebra finch chromosomes
bedtools intersect -a intermediate/freebayes/AlaArvIT_ref_AlaArv/allFreq.out -b intermediate/synteny_match/skylark_min1kb/bestMatch.list -wa -wb > intermediate/synteny_match/AlaArvIT_ref_AlaArv/allFreq.bestMatch.ZF

# Remove redundant columns
cat intermediate/synteny_match/AlaArvIT_ref_AlaArv/allFreq.bestMatch.ZF | cut -f 4,12,13,14 > intermediate/synteny_match/AlaArvIT_ref_AlaArv/allFreq.bestMatch.ZF.small

# Calculates the number of sites in 1Mbp and 100kbp windows for each sex, from that the ratio and difference in number of sites between the sexes
Rscript code/calculate_snpCount_windows.R intermediate/synteny_match/AlaArvIT_ref_AlaArv/allFreq.bestMatch.ZF.small results/AlaArvIT_ref_AlaArv/allFreq.1Mbp.out results/AlaArvIT_ref_AlaArv/allFreq.100kbp.out

