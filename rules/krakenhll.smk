KRAKENHLL = config.get('KRAKENHLL', 'krakenhll')
KRAKENHLL_BUILD = config.get('KRAKENHLL_BUILD', 'krakenhll-build')


krakenhll_all_base = expand('reports/{sample}.krakenhll.{{db}}.txt', sample=samples_all)
rule krakenhll_refseqc_all:
    input: expand(krakenhll_all_base, db='refseqc')

rule krakenhll_default_all:
    input: expand(krakenhll_all_base, db='default')

KRAKENHLL_ALL = expand(krakenhll_all_base, db=['default', 'refseqc'])
rule krakenhll_all:
    input: KRAKENHLL_ALL


KRAKENHLL_BENCHMARK_ALL = expand('benchmark/reports/{sample}.krakenhll.refseqc.txt', sample=benchmark_samples_all)
rule krakenhll_benchmark_all:
    input: KRAKENHLL_BENCHMARK_ALL

KRAKENHLL_SHELL = """
/usr/bin/time -v -o {log.time} \
{KRAKENHLL} --bzip2 --preload --db {params.db} --threads {threads}{params.options} {input} --output {params.reads} --report-file {output.report}
{PIGZ} -9 -f {params.reads}
"""

def krakenhll_db(wildcards):
    if wildcards.db == 'refseqc':
        return config['KRAKENHLL_REFSEQC_DB']
    elif wildcards.db == 'default':
        return config['KRAKENHLL_DEFAULT_DB']
    else:
        raise Exception

def krakenhll_mem(wildcards):
    if wildcards.db == 'refseqc':
        return 190
    elif wildcards.db == 'default':
        return 190
    elif wildcards.db == 'mini':
        return 4
    else:
        raise Exception

rule krakenhll:
    input: fastx_bz2_input
    output: report='reports/{seq}.krakenhll.{db}.txt',
            reads='data/{seq}.krakenhll.{db}.reads.gz'
    log: log='log/krakenhll/{seq}.{db}.log',
         time='time/krakenhll/{seq}.{db}.log'
    params: db=krakenhll_db,
            options=kraken_options,
            reads='data/{seq}.krakenhll.{db}.reads'
    threads: ALL_CORES
    resources: mem=krakenhll_mem
    shell:
        KRAKENHLL_SHELL

rule krakenhll_benchmark:
    input: fastx_bz2_input
    output: report='benchmark/reports/{seq}.krakenhll.{db}.txt',
            reads='benchmark/data/{seq}.krakenhll.{db}.reads.gz'
    log: log='benchmark/log/krakenhll/{seq}.{db}.log',
            time='benchmark/time/krakenhll/{seq}.{db}.log'
    params: db=krakenhll_db,
            options=kraken_options,
            reads='benchmark/data/{seq}.krakenhll.{db}.reads'
    threads: ALL_CORES
    resources: mem=krakenhll_mem
    benchmark: repeat('benchmark/{seq}.krakenhll.{db}.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(KRAKENHLL_SHELL, bench_record=bench_record)
        shell('truncate -s 0 {output.reads}')

rule krakenhll_refseqc_db:
    input: fna='db/refseqc/fasta/genomic.fna',
           seqmap='db/refseqc/fasta/genomic.map'
    params: dir='db/refseqc/krakenhll'
    output: 'db/refseqc/krakenhll/database.idx'
    benchmark: 'benchmark/db/krakenhll/refseqc.tsv'
    log: log='log/db/krakenhll/refseqc.log',
         time='time/db/krakenhll/refseqc.log'
    threads: ALL_CORES
    resources: mem=ALL_MEMORY
    run:
        shell('''\
        mkdir {params.dir}/taxonomy {params.dir}/library
        ln -s {TAXONOMY_DB}/names.dmp {TAXONOMY_DB}/nodes.dmp {params.dir}/taxonomy
        {PIGZ} -dc {TAXONOMY_DB}/accession2taxid/nucl_gb.accession2taxid.gz > {params.dir}/taxonomy/nucl_gb.accession2taxid
        cp {input.seqmap} {params.dir}/library/
        dropcache
        ''')

        shell('''\
        {KRAKENHLL_BUILD} --db {params.dir} --add-to-library {input.fna} > {log.log} || true  # For some reason returns 255
        /usr/bin/time -v -o {log.time} \
        {KRAKENHLL_BUILD} --db {params.dir} --build --threads {threads} --jellyfish-hash-size 20000000000 >> {log.log}
        ''', bench_record=bench_record)
