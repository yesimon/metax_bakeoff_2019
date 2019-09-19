KSLAM = config.get('KSLAM', 'SLAM')

KSLAM_DEFAULT_ALL = expand('reports/{sample}.kslam.default.tsv', sample=samples_all), expand('info/{seq}.kslam.default.classified_count.txt', seq=samples_all)
KSLAM_REFSEQC_ALL = expand('reports/{sample}.kslam.refseqc.tsv', sample=samples_all), expand('info/{seq}.kslam.refseqc.classified_count.txt', seq=samples_all)
KSLAM_ALL = KSLAM_DEFAULT_ALL + KSLAM_REFSEQC_ALL
rule kslam_all:
    input: KSLAM_ALL

rule kslam_default_all:
    input: KSLAM_DEFAULT_ALL

rule kslam_refseqc_all:
    input: KSLAM_REFSEQC_ALL

KSLAM_BENCHMARK_ALL = expand('benchmark/reports/{sample}.kslam.refseqc.tsv', sample=benchmark_samples_all)
rule kslam_benchmark_all:
    input: KSLAM_BENCHMARK_ALL

def kslam_db(wildcards):
    if wildcards.db == 'refseqc':
        return config['KSLAM_REFSEQC_DB']
    elif wildcards.db == 'default':
        return config['KSLAM_DB']
    else:
        raise Exception

KSLAM_SHELL = '''\
set +eu
source activate metax_py2
set -eu
/usr/bin/time -v -o {log.time} \
{KSLAM} --db {params.db}{params.opts} --output-file {params.out} {input} 2>&1 | tee {log.log}
cat {params.out}_PerRead | pigz -9 > {output.data}
mv {params.out}_abbreviated {output.report}
'''
rule kslam:
    input: fastq_input
    output: report='reports/{seq}.kslam.{db}.tsv',
            data='data/{seq}.kslam.{db}.tsv.gz',
    params: db=kslam_db,
            out='{seq}.out',
            opts=' --num-reads-at-once 100000000'
    log: log='log/kslam/{db}/{seq}.log',
         time='time/kslam/{db}/{seq}.log'
    shadow: 'shallow'
    threads: ALL_CORES
    resources: mem=140
    shell:
        KSLAM_SHELL

rule kslam_benchmark:
    input: fastq_input
    output: report='benchmark/reports/{seq}.kslam.{db}.tsv',
            data='benchmark/data/{seq}.kslam.{db}.tsv.gz',
    params: db=kslam_db,
            out='{seq}.out',
            opts=' --num-reads-at-once 100000000'
    log: log='benchmark/log/kslam/{db}/{seq}.log',
         time='benchmark/time/kslam/{db}/{seq}.log'
    shadow: 'shallow'
    threads: ALL_CORES
    resources: mem=140
    benchmark: repeat('benchmark/{seq}/kslam.{db}.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('{DROPCACHE}')
        shell(KSLAM_SHELL, bench_record=bench_record)

def paired_n(wildcards):
    if wildcards.seq in samples_se:
        return '1'
    elif wildcards.seq in samples_pe:
        return '2'
    else:
        return '1'

rule kslam_classified_count:
    input: 'data/{seq}.kslam.{db}.tsv.gz'
    output: 'info/{seq}.kslam.{db}.classified_count.txt'
    params: paired=paired_n
    shell:
        '''
        echo "$(zcat {input} | awk '$2 != 0' | wc -l) * {params.paired}" | bc -l > {output}
        '''

rule kslam_refseqc_db:
    input: gbff='db/refseqc/fasta/genomic.gbff'
    output: db='db/refseqc/kslam/database'
    params: tax_db=config['TAXONOMY_DB']
    log: log='log/db/kslam/refseqc.log',
         time='time/db/kslam/refseqc.log'
    benchmark: 'benchmark/db/kslam/refseqc.tsv'
    threads: ALL_CORES
    resources: mem=60
    run:
        shell('''\
        mkdir -p db/refseqc/kslam
        {DROPCACHE}
        ''')
        shell('''\
        /usr/bin/time -v -o {log.time} \
          {KSLAM} --parse-taxonomy {params.tax_db}/names.dmp {params.tax_db}/nodes.dmp --output-file taxDB | tee {log.log}
        /usr/bin/time -v -a -o {log.time} \
          {KSLAM} --output-file {output.db} --parse-genbank {input.gbff} | tee -a {log.log}
        ''', bench_record=bench_record)
