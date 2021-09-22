rule pear_merge:
    input:
        read1=get_fastq_r1,
        read2=get_fastq_r2
    output:
        assembled="pear/reads_pear_assembled.fq.gz",
        discarded="pear/reads_pear_discarded.fq.gz",
        unassembled_read1="pear/reads_pear_unassembled_r1.fq.gz",
        unassembled_read2="pear/reads_pear_unassembled_r2.fq.gz",
    log:
        'logs/pear.log'
    params:
        pval=".01",
        extra=""
    threads: threads_max
    resources:
        mem_mb=4000 # define amount of memory to be used by pear
    wrapper:
        "0.77.0/bio/pear"


