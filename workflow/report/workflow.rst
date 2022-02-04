
This HTML document describe output from the findZX pipeline. (a) This page includes settings specified in the control file, and (b) descriptions of all files listed in this document. 

- For details on methodology, see `preprint`_: 


- Instructions on how to run findZX is found on the `GitHub page`_:

- [findZX] refers to if the following snakemake file was used for the analysis: ``snakemake -s workflow/findZX``-synteny

- [findZX] refers to if the following snakemake file was used for the analysis: ``snakemake -s workflow/findZX-synteny``


----

Configuration file settings
*******************************************************

**General settings:**

- findZX analysis run name: ``{{ snakemake.config["run_name"] }}``

- Analysed samples are listed here: ``{{ snakemake.config["units"] }}``

- Reference genome: ``{{ snakemake.config["ref_genome"] }}``

- Window sizes (bp): ``{{ snakemake.config["window_sizes"] }}``

- Minimum coverage: ``{{ snakemake.config["min_cov"] }}x``

- Maximum coverage: ``{{ snakemake.config["max_cov"] }}x``

- Mismatch settings: ``{{ snakemake.config["mismatch_settings"] }}``

- Minimum size of scaffolds to plot: ``{{ snakemake.config["minSizeScaffold"] }} bp``


{% if snakemake.config["trim_reads"] %}
- Fastq files were trimmed with settings ``{{ snakemake.config["params"]["trimmomatic"]["pe"] ["trimmer"]}}``
{% elif snakemake.config["trim_and_subsample"]%}
- Fastq files were trimmed with Trimmomatic settings ``{{ snakemake.config["params"]["trimmomatic"]["pe"] ["trimmer"]}}``


- Fastq files were subsampled to ``{{ snakemake.config["subsample_basepairs"] }}`` bp
{% elif snakemake.config["subsample_only"]%}
- Fastq files were subsampled to {{ snakemake.config["subsample_only"] }}
For SNVs, the criterion ``{{ snakemake.config["params"]["trimmomatic"]["pe"] ["trimmer"]}}`` was 
{% else %}
- Fastq files were not trimmed or subsampled.
{% endif %}


**Settings for plotting [findZX]:**
{% if snakemake.config["chr_file"] %}


- List of chromosomes/scaffolds for plotting: ``({{ snakemake.config["chr_file"] }})`` 
{% endif %}
{% if snakemake.config["chr_highlight"] %}


- List of chromosomes/scaffolds to highlight (Plot type 4 only): ``({{ snakemake.config["chr_highlight"] }})`` 
{% endif %}


**Settings for plotting [findZX-synteny]:**
{% if snakemake.config["synteny_chr_file"] %}


- List of chromosomes/scaffolds for plotting: ``({{ snakemake.config["synteny_chr_file"] }})`` 
{% endif %}
{% if snakemake.config["synteny_chr_highlight"] %}

- List of chromosomes/scaffolds to highlight (Plot type 4 only): ``({{ snakemake.config["synteny_chr_highlight"] }})`` 
{% endif %}

----


Results section content
*******************************************************

This section describe the output files that are attached to this HTML report. 


**Quality control:**

{% if snakemake.config["trim_reads"] %}
`00 MultiQC reports`_: 

- MultiQC reports for trimmed and untrimmed fastq files 
{% elif snakemake.config["trim_and_subsample"]%}
`00 MultiQC reports`_: 

- MultiQC reports for trimmed and untrimmed fastq files 
{% endif %}

`01 Sample and reference genome statistics`_: 

- Reference genome statistics (calculated with `Assembly stats`_)

- Per sample heterozygosity values (total number of heterozygous sites divided by genome length)

- [findZX-synteny]: Proportion (value between 0 and 1) of 5kb windows in the study-species reference genome that was sucessfully matched with the synteny-species reference genome. 

**Output plots:**

All output plots are multi-page PDF files, where the last page contain information on what file was used to produce each plot. 

[findZX] Output plots location: ``results/{{ snakemake.config["run_name"] }}/no_synteny/plots`` 

[findZX-synteny] Output plots location:  ``results/{{ snakemake.config["run_name"] }}/synteny/{{ snakemake.config["synteny_name"] }}/plots`` 


`Output plot type 1 (genome-wide sex differences)`_: 

* Files (one for each selected window size): ``1_sexDifferences.genomeWide.{{ snakemake.config["window_sizes"] }}bp.window.pdf``

* Description: Per-sex differences in (A) heterozygosity and (B-D) genome coverage. Calculated along chromosome/scaffold positions according to the selected window sizes (bp). By default, the 50 largest scaffolds are plotted. If a list of chromosomes/scaffolds are provided in the config file ("chr_file" or "synteny_chr_file"), these are plotted.


`Output plot type 2 (genome-wide sexes separately)`_:

* Files (one for each selected window size): ``1_sexesSeparate.genomeWide.{{ snakemake.config["window_sizes"] }}bp.window.pdf``

* Description: Per-sex (A) heterozygosity and (B-D) genome coverage values (mean value for homogametic samples in purple, heterogametic samples in blue). Calculated along chromosome/scaffold positions according to the selected window sizes (bp). By default, the 50 largest scaffolds are plotted. If a list of chromosomes/scaffolds are provided in the config file ("chr_file" or "synteny_chr_file"), these are plotted. 


`Output plot type 3 (scatter plots with chromosome/scaffold length)`_:

* File name: ``3_sexDifferences.chromosome.pdf``

* Description: Scatter plots showing between-sex genome coverage and percentage of heterozygosity differences, and scaffold/chromosome length (bp). This plot is especially useful for highly fragmented reference genomes. 


`Output plot type 4 (scatter plots)`_:

* Files (one for each selected window size): ``4_sexDifferences.{{ snakemake.config["window_sizes"] }}bp.window.highlight.pdf``

* Description: Scatter plots showing between-sex genome coverage and percentage of heterozygosity differences. Dashed lines mark the genome-wide median across all genome windows. The first set of plots (page 1: A-C) show sex differences in genome coverage and heterozygosity across genome window. The second set of plots (page 2:D-F) show mean (± standard deviation) sex differences in genome coverage and heterozygosity per chromosome/scaffold, calculated from the genome windows from A-C.


`Output plot type 5 (linear models)`_:

* Files (one for each selected window size): ``5_linear_model.plot.{{ snakemake.config["window_sizes"] }}bp.pdf``

* Description: Estimates and 95% CI for linear models, testing which chromosomes siginificantly differ between sexes.


`Output plot type 6 ("confirm sexing")`_:

* Files (one for each selected mismatch setting): ``6_confirmSexing.samplesSeparately.mismatch.{{ snakemake.config["mismatch_settings"] }}.pdf`` 

* Description: These plot are based on per-individual coverage and heterozygosity values for all 5 kb windows, and can be used to (a) confirm, or identify mistakes, in the sexing of invididuals and (b) estimate suitable genome coverage filtering parameters.

----

**Output tables:**

Selected `Output tables`_ are attached to this HTML file. All output tables are found under the results directory. 

[findZX] Output tables location:  ``results/{{ snakemake.config["run_name"] }}/no_synteny/tables``

[findZX-synteny] Output tables location: ``results/{{ snakemake.config["run_name"] }}/synteny/{{ snakemake.config["synteny_name"] }}/tables``

For every selected window size, the following tables are attached:

- Test statistics from linear models, testing which chromosomes/scaffolds have sex differences (Z-transformed genome coverage and heterozygosity values) that are significantly different from zero: ``linear_model_results_estimate_CI.{{ snakemake.config["window_sizes"] }}.html``

	* Test statistics same as in `Output plot type 5 (linear models)`_

- Genome windows where genome coverage is significantly different between sexes (outside 95% CI, one file for each mismatch setting): ``diffGenomeCoverage.mismatch.{{ snakemake.config["mismatch_settings"] }}.{{ snakemake.config["window_sizes"] }}bp.outlier.out``

	* Outliers are red and blue data points in `Output plot type 1 (genome-wide sex differences)`_

- Genome windows where heterozygosity is significantly different between sexes (outside 95% CI): ``diffHeterozygosity.{{ snakemake.config["window_sizes"] }}bp.outlier.out``
 
	* Outliers are red and blue data points in `Output plot type 1 (genome-wide sex differences)`_


.. _preprint: https://www.biorxiv.org/content/10.1101/2021.10.18.464774v1
.. _GitHub page: https://github.com/hsigeman/findZX
.. _BWA mem: http://bio-bwa.sourceforge.net/
.. _Picard: https://broadinstitute.github.io/picard
.. _MultiQC: http://multiqc.info/
.. _Samtools: http://samtools.sourceforge.net/
.. _FastQC: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _Assembly stats: https://github.com/sanger-pathogens/assembly-stats