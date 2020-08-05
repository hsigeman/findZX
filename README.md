# snakemake-sex-chr

Modification of the pipline used in Sigeman et. al (2020). 

The goal of the modification is to be able to use the pipeline on more than two samples (one from each sex) and on other taxa.


# snakemake pipeline for detection of sex-linked regions
The purpose of this pipeline is to detect new and old sex-linked regions from wgs sequencing. The first version of this pipeline detected sex-linked regions for birds. Now the pipeline has been expanded to be usable on other organisms and with and without a synteny species.

The number of male and female samples does not have to be equal, the minimum requirement is one sample from each sex. 

## Requirements
Before one can run the pipeline a reference genome is needed, trimmed reads for at least one male and one female individual. Optional is a last database for a synteny species which should beplaced in a folder named data/meta/ relative to the directory from which the pipeline is ran. 

The config.txt file gives an example of how the config-file should look. Add the path to the trimmed reads, reference genome. An option is to give the path to a file with a list of chromosomes which will be used to filter out only these chromosomes in the statistical calculations. The chromsomes should be comma separated on the same line.

Create a conda environment with the environment.yml file.

The cluster.json file have to be edited if the pipeline will be ran on a cluster. Specify the account name. If a large amount of samples are used (more than 10 individuals with a genome size of 1Gbp), or an organism with a very large genome, the times and number of cores specified might have to be changed. 

### File structure
reference/genome/directory/reference.fasta  
trimmed/reads/reads.fq.gz  
data/meta/  

Make sure that the suffix of the reads and the reference genome is '.fq.gz' and '.fasta', otherwise, the pipeline will not be able to find the files. 

## How to run
The pipeline can be ran with out without a synteny species, choose the snakefile with the corresponding name (snakefile-synteny or snakefile-no-synteny).

snakemake -s snakefile-{synteny/no-synteny} -j 15 -R all --configfile config.txt --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "

-j specifies the number of jobs that can be ran simultaneously.
-R which rule to re-run, in this case it is rule all which specifies all desiered output files.

If the reference genome used is constructed from one of the individuals in the analysis, this can introduce noise from reference bias in the results, especially if the organism has a very variable genome. This can be solved by creating a consensus genome. This can be done by runing the pipeline like this:

snakemake -s snakefile-{synteny/no-synteny} -j 15 reference/genome/directory/{name_of_reference}_nonRefAf_consensus.fasta --configfile config.txt --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "

This will produce a consensus genome in the same directory as the reference genome, named the same with '_nonRefAf_consensus' added before the '.fasta' sufix. The whole pipeline can then be re-ran with the new consensus genome, change the config-file to specify to use this new reference genome and re-run the pipeline as above. The consensus genome can be created before runing the whole pipeline, or after. If it is ran before, no flagstat-files will be created since they are only specified in the rule all.

## Output
Files from all steps of the pipeline are generated.

### Figures
The result is shown in figures in the result/ folder. 
{species}.gencov_heterozygosity_indv.pdf First a figure is produced which looks at the genome coverage and heterozygosity for each individual so that one can confirm that the sexing of samples were correct.

{species}.circliz.pdf Shows difference in heterozygosity between male and female samples over chromosomes larger than 1Mbp and, if a chr.txt file were used, chromosomes in specified in the chr.txt file. Also shows the genome coverage over the same chromosomes. Three different genome coverage values are shown, 0, 2 or 4 number of mismatches allowed.

{species}.scatter.pdf A scatterplot of genome coverage vs difference in heterozygosity for 1Mbp windows colored after chromosomes. 

{species}.chr_scatter2D.pdf Scatterplots of genome coverage, difference in heterozygosity and chromsome length, where two of the variables are plotted against each other and the points are colored based on the third. Each point in these scatterplots represents a chromsome/scaffold, not a 1Mbp window as in the plots mentioned above. 

{species}.chr_scatter3D.pdf The same variables as above but in a 3D scatterplot and colored on chromosome/scaffold length.

{species}.report.html A report which states date, which samples that were used and statistics for each chromosome/scaffold.

The other files in the results directory are used to make the figures. The idea is that the figures will show if there are any sex-linked regions in ones data and if there is, one can analyse the data more troughly and make figures that suits that perticular data-set better.

### File structure
After the pipeline has ran, the following file structure will be created:

intermediate  
	bedtools  
		{reference_genome}  
			genome_5kb_windows.out  
		{species}  
			{gencov.nodup-files}  
	bwa  
		{reference_genome}  
			7bam-files}  
			{bai-files}  
			{flagstat-files}  
	freebayes  
		{reference_genome}  
			{species}.100kbp.regions  
		{species}  
			{VCF-files}  
			{difference in heterozygosity}  
			{heterozygosity for each individual}  
			{allele frequency}  
	synteny_match  
		{reference_genome}  
			genome_windows.out  
			bestMatch.list  
			bestMatch.status  
		{species}  
			{gencov.nodup.synteny-files}  
			{difference in heterozygosity.synteny-files}  
			{allele frequency.synteny-files}  
	lastal_{synteny}  
		{reference_genome}  
			{align-files}  
results  
	{species}  
		{pdf-files}  
		{1Mbp/100kbp/chr.out-files for plotting}  
		{html-report}  

# Indepth about the programs and parameters used


## Scripts in the directory code/

### calculate_hetDiff.py

### filter_allFreq.py

### heterozygosity_per_indv.py

### normalize_gencov.py and normalize_synteny_gencov.py

### make_info_html.py

### calculate_genCov_windows.R

### calculates_hetDiff_windows.R

### calculates_snpCount_windows.R

### functions.R

### histogram_indv.R

### plot_windows.R

### scatterplot_chr.R



# Solutions to known problems
In some steps a temporary directory is created. If this step fails the files are deleated but not the directory. Before this step can be re-ran one needs to delete the temporary directory, otherwise the rule will fail again.

If the sample names contain characters that can be interpetrated as an operator like '-' '+' '*' '/' the plotting_chr rule will fail as R will interpret them as an operator. This can be solved by changing the sample names.

In the filter_allele_frequency rule, the results are piped to a file with '>'. This means that the script can fail but an empty file is still created and snakemake will think that the rule ran successfully. If this happends make sure that the config file and all names are correct. The pipeline will run as it should but expects the filenames to follow these specifications. 

# Contact
Bella Sinclair, bella.sinclair@biol.lu.se


# Reference
Sigeman, H., Ponnikas, S. & Hansson, B. Whole-genome analysis across 10 songbird families within Sylvioidea reveals a novel autosome-sex chromosome fusion. Biol. Lett. 16, 20200082 (2020).

