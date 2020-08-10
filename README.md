# Snakemake Pipeline for Detection of Sex-Linked Regions
 
Modification of the pipeline used in Sigeman et. al (2020). 
 
The purpose of this pipeline is to detect new and old sex-linked regions from wgs sequencing data. The first version of this pipeline detected sex-linked regions for birds. Now the pipeline has been expanded to be usable on other organisms and with/without a synteny species.
 
The number of male and female samples does not have to be equal, the minimum requirement is one sample from each sex. 
 
 
## Requirements
Before one can run the pipeline a reference genome is needed, trimmed reads for at least one male and one female individual.  
Optional is a last database for a synteny species which should be placed in a folder named data/meta/ relative to the directory from which the pipeline is ran. 
 
The **config.txt** file gives an example of how the config-file should look. Add the path to the trimmed reads, reference genome. An option is to give the path to a file with a list of chromosomes which will be used to filter out only these chromosomes in the statistical calculations. The chromosomes should be comma separated on the same line.
 
Changed the path in the end of the **environment.yml** file to the path to your conda environment. Create a conda environment with the environment.yml file.
 
The **cluster.json** file have to be edited if the pipeline will be ran on a cluster. Specify the account name. 
If a large amount of samples are used (more than 10 individuals with a genome size of 1Gbp), or an organism with a very large genome, the times and number of cores specified might have to be changed. 
 
Make sure that the suffix of the reads and the reference genome is *'.fq.gz'* and *'.fasta'*, otherwise, the pipeline will not be able to find the files. 
 
 
## How to Run
The pipeline can be ran with and without a synteny species, choose the snakefile with the corresponding name (snakefile-synteny or snakefile-no-synteny).
 
    snakemake -s snakefile-{synteny/no-synteny} -j 15 -R all --configfile config.txt --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "
 
*-j* specifies the number of jobs that can be run simultaneously.  
*-R* specifies which rule to re-run, in this case it is rule all which specifies all desired output files.
 
If the reference genome used is constructed from one of the individuals in the analysis, this can introduce noise from reference bias in the results, especially if the organism has a very variable genome. This can be solved by creating a consensus genome. This can be done by running the pipeline like this:
 
    snakemake -s snakefile-{synteny/no-synteny} -j 15 reference/genome/directory/{name_of_reference}_nonRefAf_consensus.fasta --configfile config.txt --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "
 
This will produce a consensus genome in the same directory as the reference genome, named the same as the reference genome with *'_nonRefAf_consensus'* added before the *'.fasta'* sufix. The whole pipeline can then be re-run with the new consensus genome. Remember to change the config-file to specify this new reference genome and re-run the pipeline as above. The consensus genome can be created before running the whole pipeline, or after. If it is run before, no flagstat-files will be created since they are only specified in the rule all and here we only specify to create the consensus genome, not the files in rule all.
 
 
## Output
The result is shown in figures in the *result/* folder.   
***[species]*.gencov_heterozygosity_indv.pdf** First a figure is produced which looks at the genome coverage and heterozygosity for each individual so that one can confirm that the sexing of samples were correct.
 
***[species]*.circliz.pdf** Shows difference in heterozygosity between male and female samples over chromosomes larger than 1Mbp and, if a chr.txt file were used, chromosomes in specified in the chr.txt file. Also shows the genome coverage over the same chromosomes. Three different genome coverage values are shown, 0, 2 or 4 number of mismatches allowed.
 
***[species]*.scatter.pdf** A scatterplot of genome coverage vs difference in heterozygosity for 1Mbp windows colored after chromosomes. 
 
***[species]*.chr_scatter2D.pdf** Scatterplots of genome coverage, difference in heterozygosity and chromosome length, where two of the variables are plotted against each other and the points are colored based on the third. Each point in these scatterplots represents a chromosome/scaffold, not a 1Mbp window as in the plots mentioned above. 
 
***[species]*.chr_scatter3D.pdf** The same variables as above but in a 3D scatterplot and colored on chromosome/scaffold length.
 
***[species]*.report.html** A report which states date, which samples that were used and statistics for each chromosome/scaffold.
 
The other files in the results directory are used to make the figures. The idea is that the figures will show if there are any sex-linked regions in ones data and if there is, one can analyse the data more thoroughly and make figures that suits that particular data-set better.
 
# Contact
Bella Sinclair, bella.sinclair@biol.lu.se
 
 
# Reference
Sigeman, H., Ponnikas, S. & Hansson, B. Whole-genome analysis across 10 songbird families within Sylvioidea reveals a novel autosome-sex chromosome fusion. Biol. Lett. 16, 20200082 (2020).
