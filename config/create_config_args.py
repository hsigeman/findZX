import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--run_name", help="name of pipeline output directory (e.g. species name) - REQUIRED", required = True)
parser.add_argument("--unit_file", help="path to unit file (sample information) - REQUIRED", required = True)
parser.add_argument("--ref_genome", help="path to reference genome - REQUIRED")
parser.add_argument("--outfile", help="output config file - REQUIRED", required = True)
parser.add_argument("--threads", help="maximum number of threads (default: 2)", default="2")
parser.add_argument("--mem_max", help="maximum memory (default: 5000)", default="5000")
parser.add_argument("--syn_genome", help="path to synteny species reference genome")
parser.add_argument("--syn_species", help="synteny species name")
parser.add_argument("--window_sizes", help="window sizes (default: 100000,1000000", default="100000,1000000")
parser.add_argument("--chr_file", help="file with chromosomes to plot (optional)", default="None")
parser.add_argument("--chr_highlight", help="list chromosomes to highligt (optional). E.g. \"1, 2, X\"", default="None")
parser.add_argument("--syn_chr_file", help="file with chromosomes in synteny-species to plot (optional)", default="None")
parser.add_argument("--syn_chr_highlight", help="list synteny-species chromosomes to highligt (optional). E.g. \"1, 2, X\"", default="None")
parser.add_argument("--mismatch_settings", help="mismatch filtering settings, requires 2 arguments (default: 0.0, 0.2)", default="0.0,0.2")
parser.add_argument("--min_size_scaffold", help="minimum size scaffolds to plot (default: 10000)", default="10000")
parser.add_argument("--trim_reads", help="set to \"true\" to trim fastq files", default="false")
parser.add_argument("--trim_and_subsample", help="set to \"true\" to trim and subsample fastq files", default="false")
parser.add_argument("--subsample_only", help="set to \"true\" to subsample fastq files", default="false")
parser.add_argument("--subsample_bp", help="select number of basepairs to subsample fastq files to", default="")
parser.add_argument("--template", help="template config file (default: config/config_template.yml)", default="config/config_template.yml")

args = parser.parse_args()

dict = {"RUN_NAME": args.run_name, 
        "UNIT_FILE": args.unit_file, 
        "SYN_GENOME": args.syn_genome, 
        "THREADS": args.threads, 
        "config/dummy.fasta": args.ref_genome, 
        "SYN_SPECIES": args.syn_species,  
        "WINDOWSIZE": args.window_sizes, 
        "CHR_FILE": args.chr_file,
        "CHR_HIGHLIGHT": args.chr_highlight,
         "SYN_FILE": args.syn_chr_file,
        "SYN_HIGHLIGHT": args.syn_chr_highlight, 
        "MISMATCH": args.mismatch_settings, 
        "MINSIZESCAFFOLD": args.min_size_scaffold, 
        "TRIM_READS": args.trim_reads, 
        "TRIM_AND_SUBSAMPLE": args.trim_and_subsample, 
        "SUBSAMPLE_ONLY": args.subsample_only, 
        "SUBSAMPLE_BP": args.subsample_bp, 
        "MEM_MAX": args.mem_max}

filtered = {k: v for k, v in dict.items() if v is not None}

checkWords = list(filtered.keys())
repWords = list(filtered.values())

fin = open(args.template, "rt")
fout = open(args.outfile, "wt")

for line in fin:
	#read replace the string and write to output file
    for check, rep in zip(checkWords, repWords):
        line = line.replace(check, rep)
    fout.write(line)
#close input and output files
fin.close()
fout.close()


# python3 config/create_config_args.py --outfile hej.txt --ref_genome .test/AloPal_v1_subset.fasta --run_name monkey --unit_file .test/units.tsv --syn_genome human.fasta --syn_species homo --window_size 500000,600000