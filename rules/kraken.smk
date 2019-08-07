KRAKEN = config.get('KRAKEN', 'kraken')
KRAKEN_FILTER = config.get('KRAKEN_FILTER', 'kraken-filter')
KRAKEN_TRANSLATE = config.get('KRAKEN_TRANSLATE', 'kraken-translate')
KRAKEN_REPORT = config.get('KRAKEN_REPORT', 'kraken-report')
KRAKEN_MPA_REPORT = config.get('KRAKEN_MPA_REPORT', 'kraken-mpa-report')
KRAKEN_FILTER_THRESHOLD = config.get('KRAKEN_FILTER_THRESHOLD', 0.05)

kraken_all_base = expand('reports/{sample}.kraken.{{db}}.txt', sample=samples_all)
KRAKEN_ALL = expand(kraken_all_base, db=['default', 'refseqc'])
rule kraken_all:
    input: KRAKEN_ALL

rule kraken_default_all:
    input: expand(kraken_all_base, db='default')

rule kraken_refseqc_all:
    input: expand(kraken_all_base, db='refseqc')

KRAKEN_BENCHMARK_ALL = expand('benchmark/reports/{sample}.kraken.refseqc.txt', sample=benchmark_samples_all)
rule kraken_benchmark_all:
    input: KRAKEN_BENCHMARK_ALL

def kraken_db(wildcards):
    if wildcards.db == 'refseqc':
        return config['KRAKEN_REFSEQC_DB']
    elif wildcards.db == 'default':
        return config['KRAKEN_DEFAULT_DB']
    elif wildcards.db == 'mini':
        return config['KRAKEN_MINI_DB']
    else:
        raise Exception

def kraken_options(wildcards):
    if wildcards.seq in samples_fasta:
        return ' --fasta-input --bzip2-compressed'
    else:
        opts = ' --fastq-input --bzip2-compressed'
        if wildcards.seq in samples_pe:
            return opts + ' --paired'
        elif wildcards.seq in samples_se:
            return opts
        else:
            raise Exception

def kraken_mem(wildcards):
    if wildcards.db == 'refseqc':
        return 190
    elif wildcards.db == 'default':
        return 190
    elif wildcards.db == 'mini':
        return 4
    else:
        raise Exception


KRAKEN_SHELL = """\
echo {input}
export KRAKEN_DEFAULT_DB="{params.db}"

{{ /usr/bin/time -v -o {log.time} {KRAKEN} --preload --threads {threads}{params.options} {input} 2>&1 1>&3 3>&- | sed '/Processed [0-9]* sequences/d'; }} \
3>&1 1>&2 | tee >(pigz --best > {output.reads}) | {KRAKEN_FILTER} --threshold {KRAKEN_FILTER_THRESHOLD} | {KRAKEN_REPORT} > {output.report}
"""
rule kraken:
    input: fastx_bz2_input
    output: report='reports/{seq}.kraken.{db}.txt',
            reads='data/{seq}.kraken.{db}.reads.gz'
    log: log='log/kraken/{seq}.{db}.log',
            time='time/kraken/{seq}.{db}.log'
    params: db=kraken_db,
            options=kraken_options
    threads: ALL_CORES
    resources: mem=kraken_mem
    shell:
        KRAKEN_SHELL


rule kraken_benchmark:
    input: fastx_bz2_input
    output: report='benchmark/reports/{seq}.kraken.{db}.txt',
            reads='benchmark/data/{seq}.kraken.{db}.reads.gz'
    log: log='benchmark/log/kraken/{seq}.{db}.log',
         time='benchmark/time/kraken/{seq}.{db}.log'
    params: db=kraken_db,
            options=kraken_options
    threads: ALL_CORES
    resources: mem=kraken_mem
    benchmark: repeat('benchmark/{seq}/kraken.{db}.log', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(KRAKEN_SHELL, bench_record=bench_record)
        shell('truncate -s 0 {output}')


rule kraken_refseqc_db:
    input: fna='db/refseqc/fasta/genomic.fna'
    output: 'db/refseqc/kraken/database.idx'
    params: dir='db/refseqc/kraken',
            tax_db='db/taxonomy/20180425'
    benchmark: 'benchmark/db/kraken/refseqc.tsv'
    log: log='log/db/kraken/refseqc.log',
         time='time/db/kraken/refseqc.log'
    threads: ALL_CORES
    run:
        shell('''\
        mkdir {params.dir}/taxonomy {params.dir}/library
        ln -sf {params.tax_db}/names.dmp {params.tax_db}/nodes.dmp {params.dir}/taxonomy
        {PIGZ} -fdc {params.tax_db}/accession2taxid/nucl_gb.accession2taxid.gz > {params.dir}/taxonomy/nucl_gb.accession2taxid
        dropcache
        ''')

        shell('''\
        kraken-build --db {params.dir} --add-to-library {input.fna}
        kraken-build --db {params.dir} --build --threads {threads}
        # --jellyfish-hash-size
        ''', bench_record=bench_record)
