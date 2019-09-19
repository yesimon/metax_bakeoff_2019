METAOTHELLO = config.get('METAOTHELLO')

dbs = ['refseqc', 'default']

METAOTHELLO_ALL = expand('reports/{sample}.metaothello.{db}.tsv', sample=samples_all, db=dbs)
rule metaothello_all:
    input: METAOTHELLO_ALL

METAOTHELLO_BENCHMARK_ALL = expand('benchmark/data/{sample}.metaothello.default.tsv.gz', sample=benchmark_samples_all)
rule metaothello_benchmark_all:
    input: METAOTHELLO_BENCHMARK_ALL

def metaothello_base(db):
    if db == 'refseqc':
        base = config['METAOTHELLO_REFSEQC_DB']
    elif db == 'default':
        base = config['METAOTHELLO_DB']
    else:
        raise Exception
    return base

def metaothello_db(wildcards, input):
    return join(metaothello_base(wildcards.db), 'othello31')

def metaothello_db_map(wildcards, input):
    return join(metaothello_base(wildcards.db), 'id2tax')

def metaothello_db_names(wildcards, input):
    return join(metaothello_base(wildcards.db), 'names')

def metaothello_fx(wildcards):
    if wildcards.seq in samples_fastq:
        return 'fq'
    elif wildcards.seq in samples_fasta:
        return 'fa'
    else:
        raise Exception


METAOTHELLO_SHELL = '''\
/usr/bin/time -v -o {log.time} \
{METAOTHELLO} {params.db} {params.output_dir} 31 {threads} {params.fx} {params.paired} {params.idmap} {params.names} {input} 2>&1 | tee {log.log}
pigz --stdout -9 {params.output_dir}/taxo_assignment.txt > {output}
rm -r {params.output_dir}
'''

rule metaothello:
    input: fastx_input
    output: 'data/{seq}.metaothello.{db}.tsv.gz'
    log: log='log/metaothello/{seq}.{db}.log',
         time='time/metaothello/{seq}.{db}.log'
    params: output_dir='tmp/metaothello/{seq}',
            fx=metaothello_fx,
            paired=lambda wildcards, input: 'PE' if len(input) == 2 else 'SE',
            db=metaothello_db,
            idmap=metaothello_db_map,
            names=metaothello_db_names
    threads: ALL_CORES
    resources: mem=30
    shell:
        METAOTHELLO_SHELL

rule metaothello_benchmark:
    input: fastx_input
    output: 'benchmark/data/{seq}.metaothello.{db}.tsv.gz'
    log: log='benchmark/log/metaothello/{seq}.{db}.log',
         time='benchmark/time/metaothello/{seq}.{db}.log'
    params: output_dir='tmp/metaothello/{seq}',
            fx=metaothello_fx,
            paired=lambda wildcards, input: 'PE' if len(input) == 2 else 'SE',
            db=metaothello_db,
            idmap=metaothello_db_map,
            names=metaothello_db_names
    threads: ALL_CORES
    resources: mem=30
    benchmark: repeat('benchmark/{seq}/metaothello.{db}.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('{DROPCACHE}')
        shell(METAOTHELLO_SHELL, bench_record=bench_record)

rule metaothello_report:
    input: 'data/{seq}.metaothello.{db}.tsv.gz'
    output: report='reports/{seq}.metaothello.{db}.tsv'
    params: db='/mnt/metax/db/taxonomy/20180425'
    shell:
        '''
        metax metaothello-report --tax-dir {params.db} <({PIGZ} -dc {input}) {output}
        '''
