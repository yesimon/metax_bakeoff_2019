KRAKEN2 = config.get('KRAKEN2', 'kraken2')
KRAKEN2_FILTER_THRESHOLD = config.get('KRAKEN_FILTER_THRESHOLD', 0.05)

kraken2_all_base = expand('reports/{sample}.kraken2.{{db}}.txt', sample=samples_all)
KRAKEN2_ALL = expand(kraken2_all_base, db=['default', 'refseqc'])
rule kraken2_all:
    input: KRAKEN2_ALL

rule kraken2_default_all:
    input: expand(kraken2_all_base, db='default')

rule kraken2_refseqc_all:
    input: expand(kraken2_all_base, db='refseqc')


KRAKEN2_BENCHMARK_ALL = expand('benchmark/reports/{sample}.kraken2.refseqc.txt', sample=benchmark_samples_all)
rule kraken2_benchmark_all:
    input: KRAKEN2_BENCHMARK_ALL

def kraken2_db(wildcards):
    if wildcards.db == 'refseqc':
        return config['KRAKEN2_REFSEQC_DB']
    elif wildcards.db == 'default':
        return config['KRAKEN2_DEFAULT_DB']
    elif wildcards.db == 'mini':
        return config['KRAKEN2_MINI_DB']
    else:
        raise Exception

def kraken2_paired(wildcards):
    if wildcards.seq in samples_pe:
        return ' --paired'
    elif wildcards.seq in samples_se:
        return ''
    else:
        raise Exception

def kraken2_mem(wildcards):
    if wildcards.db == 'refseqc':
        return 190
    elif wildcards.db == 'default':
        return 190
    elif wildcards.db == 'mini':
        return 4
    else:
        raise Exception

def kraken2_options(wildcards):
    if wildcards.seq in samples_fasta:
        return ' --bzip2-compressed'
    else:
        opts = ' --bzip2-compressed'
        if wildcards.seq in samples_pe:
            return opts + ' --paired'
        elif wildcards.seq in samples_se:
            return opts
        else:
            raise Exception

KRAKEN2_SHELL = """\
echo {input}
export KRAKEN2_DEFAULT_DB="{params.db}"

/usr/bin/time -v -o {log.time} {KRAKEN2} --threads {threads}{params.options} --confidence {KRAKEN2_FILTER_THRESHOLD} --bzip2-compressed {input} --output {params.reads_tmp} --report {output.report}
{PIGZ} -f {params.reads_tmp}
"""
rule kraken2:
    input: fastq_bz2_input
    output: report='reports/{seq}.kraken2.{db}.txt',
            reads='data/{seq}.kraken2.{db}.reads.gz'
    log: log='log/kraken2/{seq}.{db}.log',
            time='time/kraken2/{seq}.{db}.log'
    params: db=kraken2_db,
            options=kraken_options,
            reads_tmp='data/{seq}.kraken2.{db}.reads'
    threads: ALL_CORES
    resources: mem=kraken2_mem
    shell:
        KRAKEN2_SHELL


rule kraken2_benchmark:
    input: fastq_bz2_input
    output: report='benchmark/reports/{seq}.kraken2.{db}.txt',
            reads='benchmark/data/{seq}.kraken2.{db}.reads.gz'
    log: log='benchmark/log/kraken2/{seq}.{db}.log',
         time='benchmark/time/kraken2/{seq}.{db}.log'
    params: db=kraken2_db,
            options=kraken_options,
            reads_tmp='benchmark/data/{seq}.kraken2.{db}.reads'
    threads: ALL_CORES
    resources: mem=kraken2_mem
    benchmark: repeat('benchmark/{seq}/kraken2.{db}.log', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(KRAKEN2_SHELL, bench_record=bench_record)
        shell('truncate -s 0 {output.reads}')


rule kraken2_refseqc_db:
    input: fna='db/refseqc/fasta/genomic.fna'
    output: 'db/refseqc/kraken2/database.idx'
    params: dir='db/refseqc/kraken2',
            tax_db='db/taxonomy/20180425'
    benchmark: 'benchmark/db/kraken2/refseqc.tsv'
    log: log='log/db/kraken2/refseqc.log',
         time='time/db/kraken2/refseqc.log'
    threads: ALL_CORES
    run:
        shell('''\
        mkdir {params.dir}/taxonomy {params.dir}/library
        ln -sf {params.tax_db}/names.dmp {params.tax_db}/nodes.dmp {params.dir}/taxonomy
        {PIGZ} -fdc {params.tax_db}/accession2taxid/nucl_gb.accession2taxid.gz > {params.dir}/taxonomy/nucl_gb.accession2taxid
        dropcache
        ''')

        shell('''\
        kraken2-build --db {params.dir} --add-to-library {input.fna}
        kraken2-build --db {params.dir} --build --threads {threads}
        ''', bench_record=bench_record)
