import sys

with open(sys.argv[1], 'r') as file:
	for line in file:
		line = line.split('.')
		sample_name = line[0]
		
		line = next(file)

		line = line.split()
		read_len = line[3]

		sys.stdout.write('%s,%s\n' % (sample_name, read_len))


