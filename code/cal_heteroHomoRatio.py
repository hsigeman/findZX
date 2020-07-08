import sys

if len(sys.argv)==1:
	sys.exit()
elif not len(sys.argv)==3:
	sys.exit()

prefix_out = sys.argv[2]

out_file_09 = open(prefix_out + '.09.bed', 'w')
out_file_08 = open(prefix_out + '.08.bed', 'w')
out_file_05 = open(prefix_out + '.05.bed', 'w')


with open(sys.argv[1], 'r') as vcfFile:
	for line in vcfFile:
		if line.startswith('#CHROM'):
			fields = line.strip('\n').split('\t')
			nr_samples_each_sex = int((len(fields) - 9)/2)
			heterogametic_sex = list(range(9, 9 + nr_samples_each_sex))
			homogametic_sex = list(range(9 + nr_samples_each_sex, 9 + 2*nr_samples_each_sex))
		elif not line.startswith('#'):
			fields = line.strip('\n').split('\t')

			heterozygout_count_hetero_sex = 0
			heterozygout_count_homo_sex = 0

			for i in heterogametic_sex:
				if fields[i].startswith('0/1'):
					heterozygout_count_hetero_sex = heterozygout_count_hetero_sex + 1
					
			heterozygot_ratio_hetero_sex = heterozygout_count_hetero_sex/nr_samples_each_sex
#			print(heterozygot_ratio_hetero_sex >= 0.9)
			if heterozygot_ratio_hetero_sex >= 0.9:
				out_file_09.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'heterogametic', str(heterozygot_ratio_hetero_sex)]) + '\n')
				out_file_08.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'heterogametic', str(heterozygot_ratio_hetero_sex)]) + '\n')
				out_file_05.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'heterogametic', str(heterozygot_ratio_hetero_sex)]) + '\n')
			elif heterozygot_ratio_hetero_sex >= 0.8:
				out_file_08.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'heterogametic', str(heterozygot_ratio_hetero_sex)]) + '\n')
				out_file_05.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'heterogametic', str(heterozygot_ratio_hetero_sex)]) + '\n')
			elif heterozygot_ratio_hetero_sex >= 0.5:
				out_file_05.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'heterogametic', str(heterozygot_ratio_hetero_sex)]) + '\n')


			for i in homogametic_sex:
				if fields[i].startswith('0/1'):
					heterozygout_count_homo_sex= heterozygout_count_homo_sex + 1

			heterozygot_ratio_homo_sex = heterozygout_count_homo_sex/nr_samples_each_sex
#			print(heterozygot_ratio_homo_sex)
			if heterozygot_ratio_homo_sex >= 0.9:
				out_file_09.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'homogametic', str(heterozygot_ratio_homo_sex)]) + '\n')
				out_file_08.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'homogametic', str(heterozygot_ratio_homo_sex)]) + '\n')
				out_file_05.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'homogametic', str(heterozygot_ratio_homo_sex)]) + '\n')
			elif heterozygot_ratio_homo_sex >= 0.8:
				out_file_08.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'homogametic', str(heterozygot_ratio_homo_sex)]) + '\n')
				out_file_05.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'homogametic', str(heterozygot_ratio_homo_sex)]) + '\n')
			elif heterozygot_ratio_homo_sex >= 0.5:
				out_file_05.write('\t'.join(fields[0:2] + [str(int(fields[1])+1), 'homogametic', str(heterozygot_ratio_homo_sex)]) + '\n')



out_file_09.close()
out_file_08.close()
out_file_05.close()



