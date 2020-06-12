import sys
import pandas as pn

if len(sys.argv)==1:
	print("Pyton-script written for the matchScaffold2Chr rule in snakemake pipe-line.")
	print("Joins two tab-separated files on the contig name and window start and stop.\n")
	sys.exit()
elif not len(sys.argv)==3:
	print("\nError:\tincorrect number of command-line arguments")
	print("Syntax:\tmatchScaffold2chr_cov.py [bestMatch.list] [genCov.out]\n")
	sys.exit()

best_match = pn.read_csv(sys.argv[1], sep='\t', header=None)
best_match = best_match.drop(best_match.columns[[3, 4, 5, 6]], axis=1)

gencov = pn.read_csv(sys.argv[2], sep='\t', header=None)

gencov_ref = pn.merge(best_match, gencov, how='inner')

sys.stdout.write(gencov_ref.to_csv(header=None, index=None, sep='\t'))

