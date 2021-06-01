if config['sra_download']:
	rule get_fastq_pe:
        output:
        	"data/{sra}_1.fastq",
        	"data/{sra}_2.fastq"
    	params:
        	# optional extra arguments
        	extra=""
    	threads: 6  # defaults to 6
    	wrapper:
        	"0.74.0/bio/sra-tools/fasterq-dump"


if config['sra_download']:
	rule get_fastq_se:
   		output:
        	"data/{sra}.fastq"
    	params:
        	extra=""
    	threads: 6
    	wrapper:
        	"0.74.0/bio/sra-tools/fasterq-dump"

