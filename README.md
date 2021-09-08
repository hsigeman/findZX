# XYZWfinder

### A snakemake pipeline for detection of sex-linked genome regions using WGS paired-end data. 



***
***


**Needed input files:**

XYZWfinder requires the following input files: 

- WGS paired-end data from at least one sample of each sex (but can handle any number of samples).

- A reference genome from the same species, generated from a homogametic individual. If no such reference genome is available, construct one from the paired-end data before starting the pipeline.

**What does XYZWfinder do?**

In essence, XYZWfinder does the following things (for details see our bioRxiv preprint at xxxxx): 

1. The paired-end data from both sexes is aligned to the reference genome

2. Genome coverage and proportion of heterozygosity per sex is calculated across 5 kb windows in the reference genome

3. The genome coverage and heterozygosity values are summarized across different (modifiable) genome window sizes, and per chromosome/scaffold. Output in the form of tables and plots. 
 
**How to run XYZWfinder**

XYZWfinder can be deployed in 2 different modes: 

- **XYZWfinder** (./workflow/XYZWfinder), in which samples are aligned to a reference genome, followed by calculation of genome coverage and heterozygosity statistics for each chromosome/scaffold as well as across genome windows of modifiable sizes. 

- **XYZWfinder-synteny** (./workflow/XYZWfinder-synteny), in which the core of the pipeline is the same as XYZWfinder but also includes a genome coordinate lift-over step between the reference genome and a second genome from another species ("synteny-species reference genome"). If the pipeline is deployed using XYZWfinder-synteny, all plots and tables will be generated based on this genome coordinate lift-over analysis. 

Details on how to install and run the pipeline is found below. 


***

# Table of Contents
1. [Installation](#installation)
2. [Verify installation by running test dataset](#test)
2. [Using the XYZWfinder workflow](#usage)
3. [Run the pipeline](#third-example)
4. [Output](#fourth-examplehttpwwwfourthexamplecom)

***    

## Installation: <a name="installation"></a>

#### **Step 1**: Obtain a copy of XYZWfinder by cloning this GitHub repository:

    git clone https://github.com/hsigeman/XYZWfinder.git
    cd XYZWfinder # Go to directory

All dependencies needed to run this workflow can be installed automatically using **conda** (see website for installation guide: https://docs.conda.io/projects/conda/en/latest/user-guide/index.html). 

#### **Step 2**: Once conda is installed, the needed software can be installed in either of two ways: 

##### Option 1: Create a minimal conda environment and install software automatically within the XYZWfinder workflow (recommended)

Enter this code to create a minimal conda environment:
 
    # The name of the environment (snakemake_basic) can be replaced with another variable
    conda create -n snakemake_basic -c conda-forge -c bioconda python=3.9.4 snakemake-wrapper-utils=0.2.0 snakemake=6.4.0 mamba=0.13.0

When launching the snakemake run, add the option: **"--use-conda"** (see [below](#test)). This enables snakemake to automatically download and install all needed software into separate conda environments for different parts of the pipeline.

##### Option 2: Install all software dependencies into a single conda environment

Use the provided conda environment file (environment.yml) to install all needed software: 

    conda env create -f environment.yml

If this option is used, omit **"--use-conda"** when launching the snakemake run.

#### **Step 3**: Once all dependencies are installed, activate the conda environment according to the instructions in the terminal. 

For example like this: 

    conda activate snakemake_basic

***

## Run the example data to make sure that all software are installed (runtime ~2 minutes per command) <a name="test"></a>

The GitHub repository contain a small test dataset (./test/Example) which can be run to verify the installation. The dataset is a small subset of paired-end WGS reads from 2 male (SRR9655170, SRR9655171) and 2 female (SRR9655168, SRR9655169) mantled howler monkeys, as well as selected scaffolds from the reference genome of this species (Genbank accession: GCA_004027835.1) and subsets of chromosomes from the human genome (Genbank accession: GCA_000001405.28).

To run the workflow using only the mantled howler monkey reference genome, run this code: 

    snakemake -s workflow/snakefile-no-synteny --configfile config/config.yml --cores 1 -R all -k --use-conda

To run the workflow with the "synteny option", run this code: 

    snakemake -s workflow/snakefile-synteny --configfile config/config.yml --cores 1 -R all -k --use-conda

*-R* specifies which rule to re-run, in this case it is rule all which specifies all desired output files.

If the workflow finishes without errors, result tables and plots will be produced here: 

    results/AloPal_test/output/

It is also possible to render an interactive HTML report using this command: 

    snakemake -s workflow/snakefile-no-synteny --configfile config/config.yml --cores 1 -R all -k --use-conda --report report.html


***

## Using the XYZWfinder workflow <a name="usage"></a>

Once the software are [installed](#installation) and [verified](#test), you can run XYZWfinder on you own dataset. To do so, carefully follow the steps outlined below. 

##### 1. Prepare input data 

i. Download WGS data from male and female individuals of the target species
ii. Download a reference genome from the expected homogametic sex of the target species (or construct one based on a WGS sample)
iii. If a "synteny-species" is used, download the reference genome of this species
iv. Create a tabular file (.tsv) with information about sample ID, heterogamety, and fastq.gz file paths (example: config/units.tsv)

##### 2. Create a configuration file

Use the template configuration files used for running the test dataset (config/config.yml) and edit where approriate. The configuration file must include the location of the tabular file containing information about the samples to be analysed (config/units.tsv; see above).



- **config/config.yml** # Specify paths to reference genome etc. 
- **config/units.tsv** # Sample information and paths to fastq files
- **config/chromosomes.list** # Optional: List of scaffolds/chromosomes in the reference genome to include in the final plots (with snakefile-no-synteny)
- **config/HS_chromosomes.list** # Optional: List of scaffolds/chromosomes in the synteny-species reference genome to include in the final plots (with snakefile-synteny)

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

Start the pipeline within a **tmux** session to ensure that the run is not stopped if you disconnect from the server (https://github.com/tmux/tmux/wiki):

    tmux new -s <name_of_session>

When using tmux to run the pipeline, make sure that the conda environment version of Python is loaded (Python 3.6) with the following command: 

    python -V # Should give: "Python 3.9.4"



If the reference genome used is constructed from one of the individuals in the analysis, this can introduce noise from reference bias in the results, especially if the organism has a very variable genome. This can be solved by creating a consensus genome. This can be done by running the pipeline like this:
 
    snakemake -s snakefile-{synteny/no-synteny} -j 15 reference/genome/directory/{name_of_reference}_nonRefAf_consensus.fasta --configfile config.txt --cluster-config cluster.json --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} "
 
This will produce a consensus genome in the same directory as the reference genome, named the same as the reference genome with *'_nonRefAf_consensus'* added before the *'.fasta'* sufix. The whole pipeline can then be re-run with the new consensus genome. Remember to change the config-file to specify this new reference genome and re-run the pipeline as above. The consensus genome can be created before running the whole pipeline, or after. If it is run before, no flagstat-files will be created since they are only specified in the rule all and here we only specify to create the consensus genome, not the files in rule all.
 
 
## Output
The result is shown in figures in the *figures/* folder.   
 

# Reference
Sigeman, H., Ponnikas, S. & Hansson, B. Whole-genome analysis across 10 songbird families within Sylvioidea reveals a novel autosome-sex chromosome fusion. Biol. Lett. 16, 20200082 (2020).

## Authors
- Hanna Sigeman (hanna.sigeman@biol.lu.se)
- Bella Sinclair (bella.sinclair@biol.lu.se)
