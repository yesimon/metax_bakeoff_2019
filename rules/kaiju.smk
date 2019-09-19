KAIJU = config.get('KAIJU', 'kaiju')
KAIJU_TABLE = config.get('KAIJU_TABLE', 'kaiju2table')

KAIJU_ALL = expand('reports/{sample}.{k}.{db}.tsv', sample=samples_all, k=['kaiju', 'kaiju_genus'], db=['nr', 'refseqc'])
rule kaiju_all:
    input: KAIJU_ALL

KAIJU_BENCHMARK_ALL = expand('benchmark/reports/{sample}.kaiju.{db}.tsv', sample=benchmark_samples_all, db='refseqc')
rule kaiju_benchmark_all:
    input: KAIJU_BENCHMARK_ALL

def kaiju_input_args(wildcards, input):
    tpl = '-i {}'.format(input[0])
    if wildcards.seq in samples_pe:
        tpl += ' -j {}'.format(input[0], input[1])
    return tpl

def kaiju_db(wildcards):
    if wildcards.db == 'nr':
        db_dir = config['KAIJU_NR_DB']
    elif wildcards.db == 'refseqc':
        db_dir = config['KAIJU_REFSEQC_DB']
    else:
        raise Exception
    return db_dir

def kaiju_db_args(wildcards):
    fmi = nodes = None
    if wildcards.db == 'nr':
        db_dir = config['KAIJU_NR_DB']
    elif wildcards.db == 'refseqc':
        db_dir = config['KAIJU_REFSEQC_DB']
    else:
        raise Exception

    for fn in os.listdir(db_dir):
        if fn.endswith('.fmi'):
            fmi = join(db_dir, fn)
        elif fn == 'nodes.dmp':
            nodes = join(db_dir, fn)
    assert fmi
    assert nodes
    return '-t {} -f {}'.format(nodes, fmi)


KAIJU_SHELL = '''\
/usr/bin/time --verbose --append -o {log.time} \
{KAIJU} {params.db_args} -m {params.min_fragment_length} {params.inf} -o {params.data} -z {threads} 2>&1 | tee {log.log}
pigz -f -9 {params.data}
'''

rule kaiju:
    input: fastx_input
    output: data='data/{seq}.kaiju.{db}.tsv.gz'
    params: db_args=kaiju_db_args,
            inf=kaiju_input_args,
            data='data/{seq}.kaiju.{db}.tsv',
            min_fragment_length=config.get('KAIJU_MIN_FRAGMENT_LENGTH', 11)
    log: log='log/kaiju/{seq}.{db}.log',
         time='time/kaiju/{seq}.{db}.log'
    resources: mem=60
    priority: 1
    threads: 4
    shell:
        KAIJU_SHELL


rule kaiju_table:
    input: 'data/{seq}.kaiju.{db}.tsv.gz'
    output: 'reports/{seq}.kaiju.{db}.tsv'
    params: kaiju_table=KAIJU_TABLE,
            db_dir=kaiju_db
    shell:
        '''
        {params.kaiju_table} -r species -e -p -v -t {params.db_dir}/nodes.dmp -n {params.db_dir}/names.dmp -o {output} <(pigz -dc {input})
        '''

rule kaiju_table_genus:
    input: 'data/{seq}.kaiju.{db}.tsv.gz'
    output: 'reports/{seq}.kaiju_genus.{db}.tsv'
    params: kaiju_table=KAIJU_TABLE,
            db_dir=kaiju_db
    shell:
        '''
        {params.kaiju_table} -r genus -e -p -v -t {params.db_dir}/nodes.dmp -n {params.db_dir}/names.dmp -o {output} <(pigz -dc {input})
        '''

KAIJU_REPORT_SHELL = '''\
{params.kaiju_table} -r species -e -p -v -t {params.db_dir}/nodes.dmp -n {params.db_dir}/names.dmp -o {output} <(pigz -dc {input})
'''
rule kaiju_benchmark:
    input: fastq_input
    output: data='benchmark/data/{seq}.kaiju.{db}.tsv.gz',
            report='benchmark/reports/{seq}.kaiju.{db}.tsv'
    params: db_dir=kaiju_db,
            db_args=kaiju_db_args,
            inf=kaiju_input_args,
            data='benchmark/data/{seq}.kaiju.{db}.tsv',
            min_fragment_length=config.get('KAIJU_MIN_FRAGMENT_LENGTH', 11)
    threads: ALL_CORES
    log: log='benchmark/log/kaiju/{seq}.{db}.log',
         time='benchmark/time/kaiju/{seq}.{db}.log'
    benchmark: repeat('benchmark/{seq}/kaiju.{db}.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('{DROPCACHE}')
        shell(KAIJU_SHELL + KAIJU_REPORT_SHELL, bench_record=bench_record)
        shell('truncate -s 0 {output}')


rule kaiju_refseqc_db:
    input: faa='db/refseqc/fasta/protein.uniq.faa.gz',
           a2t=join(TAXONOMY_DB, 'accession2taxid', 'prot.accession2taxid.vsorted')
    threads: ALL_CORES
    resources: mem=ALL_MEMORY
    shadow: 'shallow'
    benchmark: 'benchmark/db/kaiju/refseqc.tsv'
    run:
        shell(r'''\
        export LC_ALL=C

        join -t$'\t' <(zcat {input.faa} | fasta_formatter -t | sed 's/ /\t/' ) {input.a2t} | \
        awk 'BEGIN {{FS=OFS="\t"}} {{printf(">%s\n%s\n", $4, $2)}}' > protein.faa
        ''')
        shell('''\
        /mnt/metax/src/kaiju/bin/mkbwt -n {threads} -l 100000 -a protein -o proteins protein.faa
        /mnt/metax/src/kaiju/bin/mkfmi proteins
        ''', bench_record=bench_record)
        shell('rm proteins.sa proteins.bwt')
