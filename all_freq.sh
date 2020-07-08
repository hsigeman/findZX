#!/bin/bash -l
#SBATCH -A snic2019-3-374
#SBATCH -p core
#SBATCH -n 5
#SBATCH -t 1:00:00
#SBATCH -J all_freq

module load bioinfo-tools vcftools/0.1.16 BEDTools/2.29.2 python3/3.7.2

vcftools --vcf intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv.vcf --indv SJ-2333-Aarv-052-NL_S3_L001 --indv SJ-2333-Aarv-249-IT_S21_L002 --indv SJ-2333-Aarv-252-IT_S23_L002 --min-alleles 2 --max-alleles 2 --freq --remove-filtered-geno-all --minQ 20 --minDP 3 --stdout > intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-f-allFreq.bed
vcftools --vcf intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv.vcf --indv SJ-2333-Aarv-043-NL_S13_L001 --indv SJ-2333-Aarv-313-NL_S14_L001 --indv SJ-2333-Aarv-389-IT_S17_L002 --min-alleles 2 --max-alleles 2 --freq --remove-filtered-geno-all --minQ 20 --minDP 3 --stdout > intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-m-allFreq.bed

python3 code/filterAllFreq.py intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-m-allFreq.bed male > intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-allFreq.bed
python3 code/filterAllFreq.py intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-f-allFreq.bed female >> intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-allFreq.bed

sort -g -k1,1 -k2,2 -k3,3 intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-allFreq.bed > intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-allFreq-sorted.bed

bedtools intersect -a intermediate/freebayes/AlaArv_ref_AlaArv/AlaArv-allFreq-sorted.bed -b intermediate/synteny_match/AlaArv_ref_AlaArv/bestMatch.list -wa -wb > intermediate/synteny_match/AlaArv_ref_AlaArv/allFreq.bestMatch.zf

cat intermediate/synteny_match/AlaArv_ref_AlaArv/allFreq.bestMatch.zf | cut -f 4,12,13,14 > intermediate/synteny_match/AlaArv_ref_AlaArv/allFreq.bestMatch.zf.small

Rscript code/cal_snp_density_ranges.R allFreq intermediate/synteny_match/AlaArv_ref_AlaArv/allFreq.bestMatch.zf.small intermediate/synteny_match/AlaArv_ref_AlaArv/allFreq.1Mbp.out intermediate/synteny_match/AlaArv_ref_AlaArv/allFreq.100kbp.out

