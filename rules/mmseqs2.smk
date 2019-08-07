MMSEQS2 = config.get('MMSEQS2', 'mmseqs')

mmseqs2_base_all = expand('reports/{sample}.mmseqs2.{{db}}.tsv', sample=samples_all)
MMSEQS2_ALL = expand(mmseqs2_base_all, db=['nr', 'refseqc'])
rule mmseqs2_all:
    input: MMSEQS2_ALL


MMSEQS2_FASTA_ALL = expand('reports/{sample}.mmseqs2.{db}.tsv', sample=samples_fasta, db=['nr', 'refseqc'])
rule mmseqs2_fasta_all:
    input: MMSEQS2_FASTA_ALL

rule mmseqs2_refseqc_all:
    input: expand(mmseqs2_base_all, db='refseqc')

rule mmseqs2_nr_all:
    input: expand(mmseqs2_base_all, db='nr')

MMSEQS2_BENCHMARK_ALL = expand('benchmark/{sample}/mmseqs2.refseqc.tsv', sample=benchmark_samples_all)
rule mmseqs2_benchmark_all:
    input: MMSEQS2_BENCHMARK_ALL

def mmseqs2_db(wildcards):
    if wildcards.db == 'refseqc':
        return config['MMSEQS2_REFSEQC_DB']
    elif wildcards.db == 'nr':
        return config['MMSEQS2_NR_DB']
    else:
        raise Exception

def mmseqs2_tsv(wildcards):
    db = mmseqs2_db(wildcards)
    return db[:-3] + '.tsv'

MMSEQS2_SHELL = '''
/usr/bin/time -v -o {log.time} \
{MMSEQS2} createdb {input} {params.query_db} 2>&1 | tee {log.log}
/usr/bin/time -a -v -o {log.time} \
{MMSEQS2} taxonomy {params.query_db} {params.db} {params.db_tsv} {TAXONOMY_DB} {params.query_lca_db} {TMPDIR}  2>&1 | tee -a {log.log}
/usr/bin/time -a -v -o {log.time} \
{MMSEQS2} createtsv {params.query_db} {params.query_lca_db} {params.query_lca_tsv} 2>&1 | tee -a {log.log}
{PIGZ} -9 {params.query_lca_tsv}
mv {params.query_lca_tsv}.gz {output.reads}
# MMseqs doesn't automatically delete tmp workdir to cache run progress
rm -r $(readlink -f {TMPDIR}/latest) {TMPDIR}/latest
'''


MMSEQS2_NEW_SHELL = '''
/usr/bin/time -v -o {log.time} \
{MMSEQS2} createdb {input} {params.query_db} 2>&1 | tee {log.log}

/usr/bin/time -a -v -o {log.time} \
{MMSEQS2} taxonomy {params.query_db} {params.db} {params.db_tsv} {params.query_lca_db} {TMPDIR}  2>&1 | tee -a {log.log}
/usr/bin/time -a -v -o {log.time} \
{MMSEQS2} createtsv {params.query_db} {params.query_lca_db} {params.query_lca_tsv} 2>&1 | tee -a {log.log}
{MMSEQS2} taxonomyreport {params.db} {params.query_lca_db} taxonomy.report 2>&1 | tee -a {log.log}

{PIGZ} -9 {params.query_lca_tsv} > {output.reads}


# MMseqs doesn't automatically delete tmp workdir to cache run progress
rm -r $(readlink -f {TMPDIR}/latest) {TMPDIR}/latest
'''
rule mmseqs2:
    input: fastx_both_input
    output: reads='data/{seq}.mmseqs2.{db}.tsv.gz',
            report='data/{seq}.mmseqs2.{db}.report'
    log: log='log/mmseqs2/{seq}.{db}.log',
         time='time/mmseqs2/{seq}.{db}.log'
    params: db=mmseqs2_db,
            db_tsv=mmseqs2_tsv,
            query_db='queryDB',
            query_lca_db='queryLcaDB',
            query_lca_tsv='queryLca.tsv'
    shadow: 'shallow'
    threads: ALL_CORES
    resources: mem=190
    shell:
        MMSEQS2_SHELL

rule mmseqs2_benchmark:
    input: fastq_both_input
    output: reads='benchmark/data/{seq}.mmseqs2.{db}.tsv.gz'
    log: log='benchmark/log/mmseqs2/{seq}.{db}.log',
         time='benchmark/time/mmseqs2/{seq}.{db}.log'
    params: db=mmseqs2_db,
            db_tsv=mmseqs2_tsv,
            query_db='queryDB',
            query_lca_db='queryLcaDB',
            query_lca_tsv='queryLca.tsv'
    shadow: 'shallow'
    threads: ALL_CORES
    resources: mem=190
    benchmark: repeat('benchmark/{seq}/mmseqs2.{db}.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        else:
            shell('rm {params.query_lca_db} {params.query_db}')
        shell(MMSEQS2_SHELL, bench_record=bench_record)

rule mmseqs2_report:
    input: data='data/{seq}.mmseqs2.{db}.tsv.gz',
           total_reads='info/{seq}.total_reads.txt'
    output: 'reports/{seq}.mmseqs2.{db}.tsv'
    params: db=config['TAXONOMY_DB']
    shell:
        '''
        metax mmseqs-report {input.data} {output} --total-reads $(cat {input.total_reads}) --tax-dir {params.db}
        '''

rule mmseqs2_refseqc_db:
    input: faa='db/refseqc/fasta/protein.uniq.faa.gz',
           prot=join(TAXONOMY_DB, 'accession2taxid', 'prot.accession2taxid.gz'),
           dead_prot=join(TAXONOMY_DB, 'accession2taxid', 'dead_prot.accession2taxid.gz')
    output: db='db/refseqc/mmseqs2/refseqc.db'
    params: tax_db=config['TAXONOMY_DB'],
            dir='db/refseqc/mmseqs2'
    benchmark: 'benchmark/db/mmseqs2/refseqc.tsv'
    log: time='time/db/mmseqs2/refseqc.log'
    run:
      shell('''\
      /usr/bin/time -v -o {log.time} \
      {MMSEQS2} createdb {input.faa} {params.dir}/refseqc.db
      ''', bench_record=bench_record)
      # shell('''
      # metax create-mmseqs-taxonomy -o {params.dir}/refseqc.tsv --tax-db {TAXONOMY_DB} --accession2taxid {input.prot} {input.dead_prot} -- {params.dir}/refseqc.db_h
      # ''')
