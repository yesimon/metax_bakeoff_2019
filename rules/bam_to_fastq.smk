rule merge_bam_to_sorted_fastq_paired:
    input: 'bam/{sample}.bam'
    output: 'fastq/{sample}_1.fastq.bz2','fastq/{sample}_2.fastq.bz2'
    log: 'log/bam_to_fastq/{sample}.log'
    shadow: 'shallow'
    threads: ALL_CORES
    shell:
        '''
        export LC_ALL=C TMPDIR={TMPDIR}
        echo {input}

        env PICARD="{PICARD}" ~/efs/src/metax/scripts/bam_to_fastq_paired.sh {output[0]} {output[1]} {input}
        '''

def fastq_to_bam_args(wildcards, input):
    if len(input) > 1:
        return 'F1={} F2={}'.format(input[0], input[1])
    else:
        return 'F1={}'.format(input[0])

rule fastq_to_bam:
    input: fastq_input
    # output: 'tmp/bam/{seq}.bam'
    output: temp('tmp/bam/{seq}.bam')
    params: inargs=fastq_to_bam_args
    log: 'log/fastq_to_bam/{seq}.log'
    threads: ALL_CORES
    shell:
        '''
        export LC_ALL=C TMPDIR={TMPDIR}
        picard -Xmx100g -Djava.io.tmpdir={TMPDIR} FastqToSam {params.inargs} OUTPUT={output} SAMPLE_NAME={wildcards.seq} QUALITY_FORMAT=Standard MAX_RECORDS_IN_RAM=5000000
        '''
