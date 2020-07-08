#!/bin/bash -l
#SBATCH -A snic2019-3-374
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J heteroHomoRatio

# Generates files for plotting from the heterozygosity

module load bioinfo-tools vcftools/0.1.16 python3/3.7.2 BEDTools/2.29.2

# Filters VCF
vcftools --vcf ../data/AlaArv.vcf --min-alleles 2 --max-alleles 2 --remove-filtered-geno-all --minQ 20 --minDP 3 --recode --stdout > ../data/AlaArv-filtered.minQ20.minDP3.vcf

# Keeps sites if they have enought heterozygosity in one or both sexes
python3 ../code/cal_heteroHomoRatio.py ../data/AlaArv-filtered.minQ20.minDP3.vcf AlaArv-filtered.minQ20.minDP3

sort -g -k1,1 -k2,2 -k3,3 AlaArv-filtered.minQ20.minDP3.05.bed > AlaArv-filtered.minQ20.minDP3.05-sorted.bed

# Matches skylark contigs to zebra finch chromosomes
bedtools intersect -a AlaArv-filtered.minQ20.minDP3.05-sorted.bed -b ../../pipeline-test-org/intermediate/synteny_match/AlaArv_ref_AlaArv/bestMatch.list -wa -wb > AlaArv-filtered.minQ20.minDP3.05-sorted.bestMatch.zf

# Keeps certain columns
cat AlaArv-filtered.minQ20.minDP3.05-sorted.bestMatch.zf | cut -f 4,13,14,15 > AlaArv-filtered.minQ20.minDP3.05-sorted.bestMatch.zf.small

# Calculates the number of sites in 1Mbp and 100kbp windows for each sex, and from that the ratio and difference in number of sites between the sexes
Rscript code/cal_snp_density_ranges.R AlaArv-filtered.minQ20.minDP3.05 ratioHetHom/AlaArv-filtered.minQ20.minDP3.05-sorted.bestMatch.zf.small ratioHetHom/AlaArv-filtered.minQ20.minDP3.05.1Mbp.out ratioHetHom/AlaArv-filtered.minQ20.minDP3.05.100kbp.out

