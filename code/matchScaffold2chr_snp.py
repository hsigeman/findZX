import sys
import pandas as pn

if len(sys.argv)==1:
	print("Pyton-script written for the matchScaffold2Chr_snp rule in snakemake pipe-line.")
	print("Joins two tab-separated files on the contig name and window start and stop.\n")
	sys.exit()
elif not len(sys.argv)==4:
	print("\nError:\tincorrect number of command-line arguments")
	print("Syntax:\tmatchScaffold2chr_cov.py [bestMatch.list] [heterogameticSex_allDiv.out] [homogameticSex_allDiv.out]\n")
	sys.exit()

best_match = pn.read_csv(sys.argv[1], sep='\t', header=None)
best_match = best_match.drop(best_match.columns[[3, 4, 5, 6]], axis=1)

hetero_allDiv = pn.read_csv(sys.argv[2], sep='\t')
hetero_allDiv = hetero_allDiv.drop(hetero_allDiv.columns[3], axis=1)

homo_allDiv = pn.read_csv(sys.argv[3], sep='\t')
homo_allDiv = m_allDiv.drop(homo_allDiv.columns[3], axis=1)

allDiv = pn.merge(hetero_allDiv, homo_allDiv, how='inner', on=["CHROM","BIN_START","BIN_END"])

allDiv["BIN_START"] = allDiv["BIN_START"] - 1

allDiv_ref = pn.merge(best_match, allDiv, how='inner', left_on=[0,1,2], right_on=["CHROM", "BIN_START", "BIN_END"])
allDiv_ref = allDiv_ref.drop(allDiv_ref[["CHROM", "BIN_START", "BIN_END"]], axis=1)

sys.stdout.write(allDiv_ref.to_csv(header=None, index=None, sep='\t'))

