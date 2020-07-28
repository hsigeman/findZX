import sys
from datetime import date


species = argv[1]
refgenome = argv[2]
het_samples = argv[3]
homo_samples = argv[4]
nr_of_chr = argv[5]

path_statistics = argv[6]

synteny_species = argv[7]

print('<html><head><title>info_table</title></head><body><h2 style="font-family:verdana">Pipeline for detection of sex-linked regions.</h2><p style="font-family:'Courier New'">by Hanna Sigeman and Bella Sinclair, July 2020.<br>Pipeline ran:' + date.today() + '</p>')
print('<hr><p style="font-family:'Courier New'">Species:' + species + '<br>Reference genome:' + refgenome + '</p>')