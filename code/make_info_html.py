import sys
from datetime import date

if len(sys.argv) == 6:
	synteny_species = 'no'
else:
	synteny_species = sys.argv[6]

species = sys.argv[1]
refgenome = sys.argv[2]
het_samples = sys.argv[3]
homo_samples = sys.argv[4]

path_statistics = sys.argv[5]


nr_of_chr = '32'

gencov_path = 'intermediate/bedtools/' + species + '_ref_' + refgenome + '/' + species + '.gencov.nodup.nm.0.0.out'
gencov_path_synteny = 'intermediate/synteny_match/' + species + '_ref_' + refgenome + '/' + species + '.gencov.nodup.nm.0.0.' + synteny_species + '.out'


# print HTML file

print('<html>')
print('\t<head>\n\t\t<title>info_table</title>\n\t</head>')
print('\t<body>')
print('\t\t<h2 style=\"font-family:\'Arial\'\">Pipeline for detection of sex-linked regions.</h2>')
print('\t\t<p style=\"font-family:\'Arial\'\">by Hanna Sigeman and Bella Sinclair, July 2020.<br>Pipeline ran:   ' + str(date.today()) + '</p>')
print('\t\t<hr>')
print('\t\t<p style=\"font-family:\'Arial\'\">Species:   ' + species + '<br>Reference genome:   ' + refgenome + '</p>')
print('\t\t<p style=\"font-family:\'Arial\'\">Heterogametic samples:   ' + het_samples + '<br>Homogametic samples:   ' + homo_samples + '</p>')
print('\t\t<p style=\"font-family:\'Arial\'\">Synteny:   ' + synteny_species + '</p>')
print('\t\t<hr>')
print('\t\t<p style=\"font-family:\'Arial\'\">Genome coverage for each individual:</p>')
print('\t\t<p style=\"font-family:\'Courier New\'\">' + gencov_path + '</p>')

if synteny_species != 'no':
	print('\t\t<p style=\"font-family:\'Arial\'\">Genome coverage for each individual, synteny chromosomes:</p>')
	print('\t\t<p style=\"font-family:\'Courier New\'\">' + gencov_path_synteny + '</p>')

print('\t\t<p style=\"font-family:\'Arial\'\">Samples are in the order:   ' + het_samples + homo_samples + '</p>')
print('\t\t<hr>')
print('\t\t<p style=\"font-family:\'Arial\'\">Number of chromosomes/scaffolds:   ' + nr_of_chr + '</p>')
print('\t\t<hr>')
print('\t\t<h4 style=\"font-family:\'Arial\'\">Statistics calculated.</h4>')
print('\t\t<p style=\"font-family:\'Arial\'\">The ratio in genome coverage between heterogametic and homogametic individuals.<br>The difference in proportion of heterozygosity between heterogametic and homogametic individuals.</p>')

print('\t\t<table style=\"font-family:\'Courier New\'\">')

line_count = 0
with open(path_statistics, 'r') as stat_file:
	for line in stat_file:
		line_count = line_count + 1
		columns = line.strip('\n').split('\t')

		print('\t\t\t<tr>')		

		row_type = 'td'

		if line_count == 1:
			row_type = 'th'

		for col in columns:
			print('\t\t\t\t<' + row_type + '>' + col + '</' + row_type + '>')
			
		print('\t\t\t</tr>')

print('\t\t</table>')

#print('\t\t<div id="list"><p><iframe src="' + path_statistics + '" frameborder="0" height="400" width="95%"></iframe></p></div>')

print('\t\t<hr>')
print('\t\t<p style=\"font-family:\'Arial\'\">Contact: hanna.sigeman@biol.lu.se, bella.sinclair@biol.lu.se.</p>')
print('\t</body>\n</html>')



