import sys
import pandas as pn

if len(sys.argv)==1:
	sys.exit()
elif not len(sys.argv)==3:
	sys.exit()

allFreq_f = pn.read_csv(sys.argv[1], sep='\t', skiprows=1, header=None)

allFreq_f.loc[:,4]= allFreq_f.loc[:,4].str.split(':').str[1].astype(float)
allFreq_f.loc[:,5]= allFreq_f.loc[:,5].str.split(':').str[1].astype(float)

allFreq_fHet = allFreq_f.loc[(allFreq_f.loc[:,4] < 0.6) & (allFreq_f.loc[:,4] > 0.4),:]

allFreq_fHet = allFreq_fHet.drop(allFreq_fHet.columns[[2,3,4,5]], axis=1)

allFreq_fHet['pos_end'] = allFreq_f.loc[:,1]+1
allFreq_fHet['sex'] = sys.argv[2]

sys.stdout.write(allFreq_fHet.to_csv(header=None, index=None, sep='\t'))


