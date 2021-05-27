# XYZWfinder

A Snakemake pipeline for detection of sex-linked regions using WGS data. 

## Authors
- Hanna Sigeman (hanna.sigeman@biol.lu.se)
- Bella Sinclair (bella.sinclair@biol.lu.se)
 
## Usage 

### Step 1: Clone GitHub repository
    git clone https://github.com/hsigeman/XYZWfinder.git

### Step 2: Create conda environment (Tested with conda version 4.10.1)
    cd XYZWfinder
    conda env create -f environment.yml
Once all dependencies are installed, activate the conda environment according to the instructions in the terminal. 

### Step 3: Run the example data to make sure that all software are installed (runtime ~2 minutes per command)
    snakemake -s workflow/snakefile-no-synteny --configfile config/config.yml -k --cores 1 -p -R all -k
    snakemake -s workflow/snakefile-synteny --configfile config/config.yml -k --cores 1 -p -R all -k

### Step 4: Modify config files
Create configuration files with information about heterogamety of samples, and paths to files. Example config files based on a test dataset (located in .test/Example/) are here: 
- **config/config.yml** # Specify paths to reference genome etc. 
- **config/units.tsv** # Sample information and paths to fastq files
- **config/chromosomes.list** OR **config/HS_chromosomes.list** # Optional: List of chromosomes to include in the final plots
- **config/chromosomes_highlight.list** OR **config/HS_chromosomes_highlight.list** # Optional: List of chromosomes to be highlighted in the final plots

The **cluster.json** file have to be edited if the pipeline will be ran on a cluster. Specify the account name. 
If a large amount of samples are used (more than 10 individuals with a genome size of 1Gbp), or an organism with a very large genome, the times and number of cores specified might have to be changed. 
 

## Usage
The pipeline can be ran with and without a synteny species, choose the snakefile with the corresponding name (snakefile-synteny or snakefile-no-synteny).
 
    snakemake -s snakefile-{synteny/no-synteny} -j 15 -R all --configfile config.txt --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "
 
*-j* specifies the number of jobs that can be run simultaneously.  
*-R* specifies which rule to re-run, in this case it is rule all which specifies all desired output files.
 
If the reference genome used is constructed from one of the individuals in the analysis, this can introduce noise from reference bias in the results, especially if the organism has a very variable genome. This can be solved by creating a consensus genome. This can be done by running the pipeline like this:
 
    snakemake -s snakefile-{synteny/no-synteny} -j 15 reference/genome/directory/{name_of_reference}_nonRefAf_consensus.fasta --configfile config.txt --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "
 
This will produce a consensus genome in the same directory as the reference genome, named the same as the reference genome with *'_nonRefAf_consensus'* added before the *'.fasta'* sufix. The whole pipeline can then be re-run with the new consensus genome. Remember to change the config-file to specify this new reference genome and re-run the pipeline as above. The consensus genome can be created before running the whole pipeline, or after. If it is run before, no flagstat-files will be created since they are only specified in the rule all and here we only specify to create the consensus genome, not the files in rule all.
 
 
## Output
The result is shown in figures in the *figures/* folder.   
 

# Reference
Sigeman, H., Ponnikas, S. & Hansson, B. Whole-genome analysis across 10 songbird families within Sylvioidea reveals a novel autosome-sex chromosome fusion. Biol. Lett. 16, 20200082 (2020).
