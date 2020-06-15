import sys
import pandas as pn

if len(sys.argv)==1:
	print("Pyton-script written for the matchScaffold2Chr rule in snakemake pipe-line.")
	print("Joins two tab-separated files on the contig name and window start and stop.\n")
	sys.exit()
elif not len(sys.argv)==4:
	print("\nError:\tincorrect number of command-line arguments")
	print("Syntax:\tmatchScaffold2chr_cov.py [bestMatch.list] [female_allDiv.out] [male_allDiv.out]\n")
	sys.exit()

best_match = pn.read_csv(sys.argv[1], sep='\t', header=None)
best_match = best_match.drop(best_match.columns[[3, 4, 5, 6]], axis=1)

f_allDiv = pn.read_csv(sys.argv[2], sep='\t')
f_allDiv = f_allDiv.drop(f_allDiv.columns[3], axis=1)

m_allDiv = pn.read_csv(sys.argv[3], sep='\t')
m_allDiv = m_allDiv.drop(m_allDiv.columns[3], axis=1)

fm_allDiv = pn.merge(f_allDiv, m_allDiv, how='inner', on=["CHROM","BIN_START","BIN_END"])

fm_allDiv["BIN_START"] = fm_allDiv["BIN_START"] - 1

fm_allDiv_ref = pn.merge(best_match, fm_allDiv, how='inner', left_on=[0,1,2], right_on=["CHROM", "BIN_START", "BIN_END"])
fm_allDiv_ref = fm_allDiv_ref.drop(fm_allDiv_ref[["CHROM", "BIN_START", "BIN_END"]], axis=1)

sys.stdout.write(fm_allDiv_ref.to_csv(header=None, index=None, sep='\t'))

