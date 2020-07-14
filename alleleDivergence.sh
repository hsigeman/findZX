#!/bin/bash -l
#SBATCH -A snic2019-3-374
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J all_div

# Generates a file from which the allele divergence can be plotted.


module load bioinfo-tools vcftools/0.1.16 BEDTools/2.29.2 python3/3.7.2

# Allele divergence for heterogametic and homogametic samples
#vcftools --vcf intermediate/freebayes/AlaArvIT_ref_AlaArv/AlaArvIT.biallelic.minQ20.minDP3.vcf --window-pi 5000 --indv SJ-2333-Aarv-249-IT_S21_L002 --indv SJ-2333-Aarv-441-IT_S25_L002 --indv SJ-2333-Aarv-460-IT_S27_L002 --stdout > intermediate/freebayes/AlaArvIT_ref_AlaArv/allDiv.heterogametic.out
#vcftools --vcf intermediate/freebayes/AlaArvIT_ref_AlaArv/AlaArvIT.biallelic.minQ20.minDP3.vcf --window-pi 5000 --indv SJ-2333-Aarv-389-IT_S17_L002 --indv SJ-2333-Aarv-449-IT_S19_L002 --indv SJ-2333-Aarv-939-IT_S20_L002 --stdout > intermediate/freebayes/AlaArvIT_ref_AlaArv/allDiv.homogametic.out

# Matches skylark scaffolds too zebra finch chromosomes
python3 code/matchScaffold2chr_snp.py intermediate/synteny_match/skylark_min1kb/bestMatch.list intermediate/freebayes/AlaArvIT_ref_AlaArv/allDiv.heterogametic.out intermediate/freebayes/AlaArvIT_ref_AlaArv/allDiv.homogametic.out > intermediate/synteny_match/AlaArvIT_ref_AlaArv/allDiv.bestMatch.ZF.out

# Calculates the pi ratio between the sexes
Rscript code/calculate_gencov_windows.R intermediate/synteny_match/AlaArvIT_ref_AlaArv/allDiv.bestMatch.ZF.out results/AlaArvIT_ref_AlaArv/AlaArvIT.allDiv.1Mbp.out results/AlaArvIT_ref_AlaArv/AlaArvIT.allDiv.100kbp.out

