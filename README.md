# snakemake-sex-chr

Modification of the pipline used in Sigeman et. al (2020). 

The goal of the modification is to be able to use the pipeline on more than two samples (one from each sex) and on other taxa.

Contact: Bella Sinclair, bella.sinclair@biol.lu.se

Reference:
Sigeman, H., Ponnikas, S. & Hansson, B. Whole-genome analysis across 10 songbird families within Sylvioidea reveals a novel autosome-sex chromosome fusion. Biol. Lett. 16, 20200082 (2020).



Work directory
	- intermediate
		- bedtools
			- species_refSpecies
				- gencov-files
		- bwa
			- species_refSpecies
				- bam-files (nodup, nm)
				- bai-files (nodup, nm)
				- flagstat-files
				- status-file
				- duplicatedata.txt
		- freebayes
			- species_refSpecies
				- VCF-file
				- hetero/homogametic.bed-files (allele divergence, allele frequency, proportion heterozygosity)
				- status-file
				- regions-file
		- lastal_ZF
			- align-files
		- synteny_match
			- species_refSpecies
				- bestmatch-filer
				- gencov.zf-filer
				- allDiv, allFreq, propHetZyg
	- results
		- species_refSpecies
			- plots
			- 100kbp-filer (gencov, allDiv, allFreq, propHetZyg)
			- 1Mbp-filer (gencov, allDiv, allFreq, propHetZyg)
	- code
	- data
		- genome
			- consensus genome
		- meta
			- zf-database
Reads directory
	- fastQ-files
Reference genome
	- reference genome
	- index files

