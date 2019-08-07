METAPHLAN2 = config.get('METAPHLAN2', 'metaphlan2.py')
MPA_DIR = config.get('MPA_DIR')

METAPHLAN2_ALL = expand('reports/{sample}.metaphlan2.tsv', sample=samples_all) + expand('classified_count/{sample}.metaphlan2.txt', sample=samples_all)
rule metaphlan2_all:
    input: METAPHLAN2_ALL

METAPHLAN2_BENCHMARK_ALL = expand('benchmark/reports/{sample}.metaphlan2.tsv', sample=benchmark_samples_all)
rule metaphlan2_benchmark_all:
    input: METAPHLAN2_BENCHMARK_ALL

METAPHLAN2_SHELL = '''\
# Don't use bowtie caching functionality because it's so fast anyways
rm -f {output.bowtie2out}

source activate metax_py2
/usr/bin/time -v -o {log.time} \
{METAPHLAN2} --input_type fastq --nproc {threads} \
--bowtie2out {output.bowtie2out} \
--mpa_pkl {MPA_DIR}/mpa_v20_m200.pkl --bowtie2db {MPA_DIR}/mpa_v20_m200 \
{input} {output.report} 2>&1 | tee {log.log}
'''
rule metaphlan2:
    input: fastq_both_input
    output: report='reports/{seq}.metaphlan2.tsv',
            bowtie2out='data/{seq}.bowtie2out.txt'
    log: log='log/metaphlan2/{seq}.log',
         time='time/metaphlan2/{seq}.log'
    resources: mem=4
    threads: 1
    shadow: 'shallow'
    shell:
        METAPHLAN2_SHELL


rule metaphlan2_benchmark:
    input: fastq_both_input
    output: report='benchmark/reports/{seq}.metaphlan2.tsv',
            bowtie2out='benchmark/data/{seq}.bowtie2out.txt'
    log: log='benchmark/log/metaphlan2/{seq}.log',
         time='benchmark/time/metaphlan2/{seq}.log'
    resources: mem=4
    threads: ALL_CORES
    shadow: 'shallow'
    benchmark: repeat('benchmark/{seq}/metaphlan2.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(METAPHLAN2_SHELL, bench_record=bench_record)


rule metaphlan2_classified_count:
    input: 'data/{seq}.bowtie2out.txt'
    output: 'classified_count/{seq}.metaphlan2.txt'
    shell:
        '''
        cut -f1 {input} | uniq | wc -l > {output}
        '''
