import sys

if len(sys.argv)==1:
	print("\nFilters an allele frequency file and keeps sites that have a")
	print("frequency between 0.4 and 0.6. Prints to stdout.\n")
	sys.exit()
elif not len(sys.argv)==3:
	print("\nERROR: wrong number of input arguments")
	print("Call: python3 {script} {Allele Frequency} {sex}\n")
	sys.exit()


sex = sys.argv[2]

with open(sys.argv[1], 'r') as allFreq:

	for line in allFreq:

		if not line.startswith('CHROM\tPOS'):

			fields = line.split('\t')
			frequency = float(fields[4].split(':')[1])

			if frequency < 0.6 and frequency > 0.4:

				sys.stdout.write('\t'.join([fields[0], fields[1], str(int(fields[1]) + 1), sex]) + '\n')


