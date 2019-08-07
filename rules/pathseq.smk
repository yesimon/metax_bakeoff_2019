from os.path import join
GATK = config.get('GATK')
PATHSEQ_DB = config.get('PATHSEQ_DB')

PATHSEQ_ALL = expand('reports/{seq}.pathseq.txt', seq=samples_all), expand('classified_count/{seq}.pathseq.txt', seq=samples_all)
rule pathseq_all:
    input: PATHSEQ_ALL

PATHSEQ_BENCHMARK_ALL = expand('benchmark/reports/{seq}.pathseq.txt', seq=benchmark_samples_all)
rule pathseq_benchmark_all:
    input: PATHSEQ_BENCHMARK_ALL

PATHSEQ_SHELL = '''\
/usr/bin/time -v -o {log.time} \
gatk --java-options -Djava.io.tmpdir={TMPDIR} PathSeqPipelineSpark --input {input} --microbe-bwa-image {params.bwa_db} --microbe-fasta {params.microbe_fasta} --taxonomy-file {params.taxonomy} --scores-output {output.report} --output {output.reads} 2>&1 | tee {log.log}
'''
rule pathseq:
    input: 'tmp/bam/{seq}.bam'
    output: report='reports/{seq}.pathseq.txt',
            reads='data/{seq}.pathseq.bam'
    log: log='log/pathseq/{seq}.log',
         time='time/pathseq/{seq}.log'
    params: bwa_db=join(PATHSEQ_DB, 'pathseq_microbe.fa.img'),
            microbe_fasta=join(PATHSEQ_DB, 'pathseq_microbe.fa'),
            taxonomy=join(PATHSEQ_DB, 'pathseq_taxonomy.db'),
    threads: ALL_CORES
    resources: mem=170
    shell:
        PATHSEQ_SHELL

rule pathseq_benchmark:
    input: 'tmp/bam/{seq}.bam'
    output: report='benchmark/reports/{seq}.pathseq.txt',
            reads='benchmark/data/{seq}.pathseq.bam'
    log: log='benchmark/log/pathseq/{seq}.log',
         time='benchmark/time/pathseq/{seq}.log'
    params: bwa_db=join(PATHSEQ_DB, 'pathseq_microbe.fa.img'),
            microbe_fasta=join(PATHSEQ_DB, 'pathseq_microbe.fa'),
            taxonomy=join(PATHSEQ_DB, 'pathseq_taxonomy.db'),
    threads: ALL_CORES
    resources: mem=170
    benchmark: repeat('benchmark/{seq}/pathseq.log', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(PATHSEQ_SHELL, bench_record=bench_record)

def is_paired(wildcards):
    if wildcards.seq in samples_pe:
        return True
    else:
        return False

rule pathseq_classified_count:
    input: 'data/{seq}.pathseq.bam'
    output: 'classified_count/{seq}.pathseq.txt'
    params: is_paired=is_paired
    run:
        if params.is_paired:
            shell('expr "$(samtools view -c -f 0x41 -F 0xf04 {input})" + "$(samtools view -c -f 0x81 -F 0xf04 {input})" > {output}')
        else:
            shell('samtools view -c -F 0xf04 {input} > {output}')
