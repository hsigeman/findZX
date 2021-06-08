# XYZWfinder

A Snakemake pipeline for detection of sex-linked regions using WGS data. 

## Authors
- Hanna Sigeman (hanna.sigeman@biol.lu.se)
- Bella Sinclair (bella.sinclair@biol.lu.se)
 
## Install software and run XYZWfinder with example dataset:

### Step 1: Clone GitHub repository
    git clone https://github.com/hsigeman/XYZWfinder.git
    cd XYZWfinder # Go to directory

### Step 2: Create conda environment (Tested with conda version 4.10.1)

We recommend the user to create a minimal conda environment using the following code, and to launch the snakemake run with the following additional options: **"--use-conda --conda-frontend mamba"** (see Step 3). This downloads and installs separate conda environments for different parts of the pipeline.

    conda create -n snakemake_basic -c conda-forge python=3.9.4 snakemake-wrapper-utils=0.2.0 snakemake=6.4.0 mamba=0.13.0

Alternatively, all software dependencies can be installed using the provided conda environment file, in which case **"--use-conda --conda-frontend mamba"** should be omitted when launching the snakemake run: 

    conda env create -f environment.yml

Once all dependencies are installed, activate the conda environment according to the instructions in the terminal. 

### Step 3: Run the example data to make sure that all software are installed (runtime ~2 minutes per command)
    snakemake -s workflow/snakefile-no-synteny --configfile config/config.yml --cores 1 -R all -k --use-conda --conda-frontend mamba
    snakemake -s workflow/snakefile-synteny --configfile config/config.yml --cores 1 -R all -k --use-conda --conda-frontend mamba

## Create configuration files: 

To run XYZWfinder on your own dataset, you need to modify or create new configuration files. Use the template configuration files used for running the test dataset (config/config.yml) and edit where approriate. The configuration file must include the location of a tabular file containing information about the samples to be analysed (for the example dataset: config/units.tsv).


**Sample information (config/units.tsv):** 
| sample            | unit          | fq1                                     | fq2                                     |
|-------------------|---------------|-----------------------------------------|-----------------------------------------|
| subset_SRR9655168 | homogametic   | .test/Example/subset_SRR9655168_1.fq.gz | .test/Example/subset_SRR9655168_2.fq.gz |
| subset_SRR9655169 | homogametic   | .test/Example/subset_SRR9655169_1.fq.gz | .test/Example/subset_SRR9655169_2.fq.gz |
| subset_SRR9655170 | heterogametic | .test/Example/subset_SRR9655170_1.fq.gz | .test/Example/subset_SRR9655170_2.fq.gz |
| subset_SRR9655171 | heterogametic | .test/Example/subset_SRR9655171_1.fq.gz | .test/Example/subset_SRR9655171_2.fq.gz |



### Step 1: Modify config files
Create configuration files with information about heterogamety of samples, and paths to files. Example config files based on a test dataset (located in .test/Example/) are here: 
- **config/config.yml** # Specify paths to reference genome etc. 
- **config/units.tsv** # Sample information and paths to fastq files
- **config/chromosomes.list** # Optional: List of scaffolds/chromosomes in the reference genome to include in the final plots (with snakefile-no-synteny)
- **config/HS_chromosomes.list** # Optional: List of scaffolds/chromosomes in the synteny-species reference genome to include in the final plots (with snakefile-synteny)
- **config/chromosomes_highlight.list** # Optional: List of scaffolds/chromosomes in the reference genome to be highlighted in the final plots (with snakefile-no-synteny)
- **config/HS_chromosomes_highlight.list** # Optional: List of scaffolds/chromosomes in the synteny-species reference genome to be highlighted in the final plots (with snakefile-synteny)

The **cluster.json** file have to be edited if the pipeline will be ran on a cluster. Specify the account name. 
If a large amount of samples are used (more than 10 individuals with a genome size of 1Gbp), or an organism with a very large genome, the times and number of cores specified might have to be changed. 
 

## Run the pipeline:
The pipeline can be run with and without a synteny species. Choose the snakefile with the corresponding name (snakefile-synteny or snakefile-no-synteny).
    
The first step of the pipeline is optional trimming of all samples, with fastqc and multiqc being run on the samples before and after trimming. To only run this part of the pipeline (if TRIM_SAMPLES in the config file is set to "TRUE"), run the pipeline like this: 

    snakemake -s workflow/snakefile-{synteny/no-synteny} --configfile config/config.yml -k --cores {N} --use-conda --conda-frontend mamba -R multiqc_stop --notemp

If the multiqc report shows that the trimming was successful (or if trimming is not needed), run the pipeline again with this command:

    snakemake -s workflow/snakefile-{synteny/no-synteny} --configfile config/config.yml -k --cores {N} --use-conda --conda-frontend mamba -R all


If the pipeline is run on a server cluster (e.g. SLURM), a configuration file is needed (example cluster.yaml), and the command to start the pipeline should be written like this: 

    snakemake -s workflow/snakefile-{synteny/no-synteny} -j 15 -R all --configfile config/config.yml --cluster-config cluster.yaml --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "
 
*-j* specifies the number of jobs that can be run simultaneously.  
*-R* specifies which rule to re-run, in this case it is rule all which specifies all desired output files.
 





If the reference genome used is constructed from one of the individuals in the analysis, this can introduce noise from reference bias in the results, especially if the organism has a very variable genome. This can be solved by creating a consensus genome. This can be done by running the pipeline like this:
 
    snakemake -s snakefile-{synteny/no-synteny} -j 15 reference/genome/directory/{name_of_reference}_nonRefAf_consensus.fasta --configfile config.txt --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "
 
This will produce a consensus genome in the same directory as the reference genome, named the same as the reference genome with *'_nonRefAf_consensus'* added before the *'.fasta'* sufix. The whole pipeline can then be re-run with the new consensus genome. Remember to change the config-file to specify this new reference genome and re-run the pipeline as above. The consensus genome can be created before running the whole pipeline, or after. If it is run before, no flagstat-files will be created since they are only specified in the rule all and here we only specify to create the consensus genome, not the files in rule all.
 
 
## Output
The result is shown in figures in the *figures/* folder.   
 

# Reference
Sigeman, H., Ponnikas, S. & Hansson, B. Whole-genome analysis across 10 songbird families within Sylvioidea reveals a novel autosome-sex chromosome fusion. Biol. Lett. 16, 20200082 (2020).
