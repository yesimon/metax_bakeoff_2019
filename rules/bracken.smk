BRACKEN = config.get('BRACKEN', 'est_abundance.py')

def bracken_db(wildcards):
    if wildcards.db == 'default':
        return config.get('BRACKEN_DEFAULT_DB')
    elif wildcards.db == 'refseqc':
        return config.get('BRACKEN_REFSEQC_DB')
    else:
        raise Exception

BRACKEN_SHELL = '''\
set +eu
source activate metax_py2
set -eu
/usr/bin/time -v -o {log.time} \
{BRACKEN} -i {input} -k {params.db} -o {output} 2>&1 | tee {log.log}
'''

rule bracken:
    input: 'reports/{seq}.kraken.{db}.txt'
    output: 'reports/{seq}.bracken.{db}.txt'
    params: db=bracken_db
    log: log='log/bracken/{seq}.{db}.txt',
         time='time/bracken/{seq}.{db}.txt'
    shell:
        BRACKEN_SHELL

rule bracken_genus:
    input: 'reports/{seq}.kraken.{db}.txt'
    output: 'reports/{seq}.bracken_genus.{db}.txt'
    params: exe=BRACKEN,
            db=bracken_db,
            new_report='reports/{seq}.kraken.{db}_bracken.txt'
    log: log='log/bracken/genus/{seq}.{db}.txt',
         time='time/bracken/genus/{seq}.{db}.txt'
    shell:
        '''
        export PS1=
        source activate metax_py2
        echo {input}
        /usr/bin/time -v -o {log.time} \
        {params.exe} -i {input} -k {params.db} -o {output} -l G 2>&1 | tee {log.log}
        '''

BRACKEN_ALL = expand('reports/{sample}.bracken{rank}.{db}.txt', sample=samples_all, db=['default', 'refseqc'], rank=['', '_genus'])
rule bracken_all:
    input: BRACKEN_ALL

BRACKEN_BENCHMARK_ALL = expand('benchmark/reports/{sample}.bracken.refseqc.txt', sample=benchmark_samples_all)

rule bracken_benchmark_all:
    input: BRACKEN_BENCHMARK_ALL

rule bracken_benchmark:
    input: 'benchmark/reports/{seq}.kraken.refseqc.txt'
    output: 'benchmark/reports/{seq}.bracken.refseqc.txt'
    params: db='db/refseqc/bracken/KMER_DISTR.150.txt'
    log: log='benchmark/log/bracken/{seq}.refseqc.txt',
        time='benchmark/time/bracken/{seq}.refseqc.txt'
    benchmark: repeat('benchmark/{seq}/bracken.refseqc.tsv', 2)
    threads: ALL_CORES
    run:
        if benchmark_i == 0:
            shell('{DROPCACHE}')
        shell(BRACKEN_SHELL, bench_record=bench_record)

rule bracken_refseqc_db:
    input: 'db/refseqc/kraken/database.kdb'
    output: 'db/refseqc/bracken/KMER_DISTR.150.txt'
    params: kraken_dir='db/refseqc/kraken',
            dir='db/refseqc/bracken',
            tax_db='db/taxonomy/20180425'
    benchmark: 'benchmark/db/bracken/refseqc.tsv'
    log: log='log/db/bracken/refseqc.log',
        time='time/db/bracken/refseqc.log'
    threads: ALL_CORES
    run:
        shell('''\
        {DROPCACHE}
        ''')
        shell(r'''\
        export PS1=
        source activate metax_py2
        /usr/bin/time -v -o {log.time} \
        kraken --db {params.kraken_dir} --fasta-input --threads 32 <( find -L {params.kraken_dir}/library \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) -exec cat {{}} + )  > {params.dir}/database.kraken
        /usr/bin/time -v -a -o {log.time} \
        perl ~/efs/metax/Bracken/src/count-kmer-abundances.pl --db {params.kraken_dir} --read-length 150 {params.dir}/database.kraken > {params.dir}/database150mers.kraken_cnts
        /usr/bin/time -v -a -o {log.time} \
        python ~/efs/metax/Bracken/src/generate_kmer_distribution.py -i {params.dir}/database150mers.kraken_cnts -o {params.dir}/KMER_DISTR.150.txt
        ''', bench_record=bench_record)
