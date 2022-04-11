import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument("-f_hom", "--forward_homogametic", help="REQUIRED", required = True)
parser.add_argument("-r_hom", "--reverse_homogametic", help="REQUIRED", required = True)
parser.add_argument("-f_het", "--forward_heterogametic", help="REQUIRED", required = True)
parser.add_argument("-r_het", "--reverse_heterogametic", help="REQUIRED", required = True)
parser.add_argument("-f_suffix", "--forward_suffix", help="REQUIRED. Suffix for SAMPLE_1.fq.gz is \"_1.fq.gz\"", required = True)
parser.add_argument("-out", "--output_file", help="REQUIRED. Output file", required = True)

args = parser.parse_args()


forward_homogametic = pd.read_table(args.forward_homogametic, header=None)
forward_homogametic.columns = ['fq1']
forward_homogametic['group'] = 'monogynous'

reverse_homogametic = pd.read_csv(args.reverse_homogametic, header=None)
reverse_homogametic.columns = ['fq2']

forward_heterogametic = pd.read_csv(args.forward_heterogametic, header=None)
forward_heterogametic.columns = ['fq1']
forward_heterogametic['group'] = "polygynous"

reverse_heterogametic = pd.read_csv(args.reverse_heterogametic, header=None)
reverse_heterogametic.columns = ['fq2']


forward_homogametic['sample'] = forward_homogametic["fq1"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
forward_homogametic['sample'] = forward_homogametic['sample'].str.replace(args.forward_suffix, "", regex = True)
forward_homogametic = forward_homogametic[["sample", "group", "fq1"]]
 
forward_heterogametic['sample'] = forward_heterogametic["fq1"].apply(lambda x: os.path.splitext(os.path.basename(x))[0])
forward_heterogametic['sample'] = forward_heterogametic['sample'].str.replace(args.forward_suffix, "", regex=True)
forward_heterogametic = forward_heterogametic[["sample", "group", "fq1"]]


output1 = pd.concat([forward_homogametic, reverse_homogametic], axis=1)
output2 = pd.concat([forward_heterogametic, reverse_heterogametic], axis=1)
output = pd.concat([output1, output2], axis = 0)
output = pd.DataFrame(output) 

#output = output.to_string(index=False)
output.to_csv(args.output_file, sep='\t', index=False, header=True)
#print(output)