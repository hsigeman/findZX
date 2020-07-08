#!/bin/bash -l
#SBATCH -A snic2019-3-374
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J heteroHomoRatio

module load bioinfo-tools vcftools/0.1.16 python3/3.7.2 BEDTools/2.29.2

vcftools --vcf ../data/AlaArv.vcf --min-alleles 2 --max-alleles 2 --remove-filtered-geno-all --minQ 20 --minDP 3 --recode --stdout > ../data/AlaArv-filtered.minQ20.minDP3.vcf

python3 ../code/cal_heteroHomoRatio.py ../data/AlaArv-filtered.minQ20.minDP3.vcf AlaArv-filtered.minQ20.minDP3

sort -g -k1,1 -k2,2 -k3,3 AlaArv-filtered.minQ20.minDP3.05.bed > AlaArv-filtered.minQ20.minDP3.05-sorted.bed

bedtools intersect -a AlaArv-filtered.minQ20.minDP3.05-sorted.bed -b ../../pipeline-test-org/intermediate/synteny_match/AlaArv_ref_AlaArv/bestMatch.list -wa -wb > AlaArv-filtered.minQ20.minDP3.05-sorted.bestMatch.zf

cat AlaArv-filtered.minQ20.minDP3.05-sorted.bestMatch.zf | cut -f 4,13,14,15 > AlaArv-filtered.minQ20.minDP3.05-sorted.bestMatch.zf.small

Rscript code/cal_snp_density_ranges.R AlaArv-filtered.minQ20.minDP3.05 ratioHetHom/AlaArv-filtered.minQ20.minDP3.05-sorted.bestMatch.zf.small ratioHetHom/AlaArv-filtered.minQ20.minDP3.05.1Mbp.out ratioHetHom/AlaArv-filtered.minQ20.minDP3.05.100kbp.out

