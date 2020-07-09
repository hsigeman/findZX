import sys
import pandas as pn

if len(sys.argv)==1:
	print("\nFilters an allele frequency file and keeps sites that have a")
	print("frequency between 0.4 and 0.6. Prints to stdout.\n")
	sys.exit()
elif not len(sys.argv)==3:
	print("\nERROR: wrong number of input arguments")
	print("Call: python3 {script} {Allele Frequency} {sex}\n")
	sys.exit()

allFreq = pn.read_csv(sys.argv[1], sep='\t', skiprows=1, header=None)

allFreq.loc[:,4]= allFreq.loc[:,4].str.split(':').str[1].astype(float)

allFreqHet = allFreq.loc[(allFreq.loc[:,4] < 0.6) & (allFreq.loc[:,4] > 0.4),:]

allFreqHet = allFreqHet.drop(allFreqHet.columns[[2,3,4,5]], axis=1)

allFreqHet['pos_end'] = allFreq.loc[:,1]+1
allFreqHet['sex'] = sys.argv[2]

sys.stdout.write(allFreqHet.to_csv(header=None, index=None, sep='\t'))


