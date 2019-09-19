TAXMAPS = config.get('TAXMAPS', 'taxMaps')
GEMPATH = config.get('GEMPATH')

TAXMAPS_DEFAULT_ALL = expand('reports/{sample}.taxmaps.default.tsv', sample=samples_se + samples_pe)
TAXMAPS_REFSEQC_ALL = expand('reports/{sample}.taxmaps.refseqc.tsv', sample=samples_se + samples_pe)

TAXMAPS_ALL = TAXMAPS_DEFAULT_ALL + TAXMAPS_REFSEQC_ALL
rule taxmaps_default_all:
    input: TAXMAPS_DEFAULT_ALL
rule taxmaps_refseqc_all:
    input: TAXMAPS_REFSEQC_ALL
rule taxmaps_all:
    input: TAXMAPS_ALL

TAXMAPS_BENCHMARK_ALL = expand('benchmark/reports/{sample}.taxmaps.refseqc.tsv', sample=benchmark_samples_se + benchmark_samples_pe)
rule taxmaps_benchmark_all:
    input: TAXMAPS_BENCHMARK_ALL

def taxmaps_db(wildcards):
    if wildcards.db == 'default':
        return config.get('TAXMAPS_DB')
    elif wildcards.db == 'refseqc':
        return config.get('TAXMAPS_REFSEQC_DB')
    else:
        raise Exception

def taxmaps_tax(wildcards):
    if wildcards.db == 'default':
        return config.get('TAXMAPS_TAX')
    elif wildcards.db == 'refseqc':
        return config.get('TAXMAPS_REFSEQC_TAX')
    else:
        raise Exception

taxmaps_rule = '''\
set +eu
source activate metax_py2
set -eu
rm -rf tmp/taxmaps_{wildcards.seq}
/usr/bin/time -v -o {log.time} \
  env PATH={GEMPATH}:$PATH \
  {TAXMAPS} -f {input} -d {params.db} -c {threads} -t {params.tax} -p sample -o tmp/taxmaps_{wildcards.seq}
  {PIGZ} -9 -c tmp/taxmaps_{wildcards.seq}/txM.sample.map/sample.merged.map.lca > {output.lca}
  mv tmp/taxmaps_{wildcards.seq}/txM.sample.out/sample.merged.map.lca.summary {output.summary}
  mv tmp/taxmaps_{wildcards.seq}/txM.sample.out/sample.tbl {output.report}
  mv tmp/taxmaps_{wildcards.seq}/txM.sample.log/sample.prins.err {log.prins}
  mv tmp/taxmaps_{wildcards.seq}/txM.sample.log/sample.*.gem.err {log.gem}
  rm -r tmp/taxmaps_{wildcards.seq}
'''

rule taxmaps:
    input: fastq_both_input
    output: report='reports/{seq}.taxmaps.{db}.tsv',
            summary='reports/extra/{seq}.taxmaps.lca.summary.{db}.tsv',
            lca='data/{seq}.taxmaps.lca.{db}.tsv'
    priority: 2
    log: prins='log/taxmaps/{seq}.{db}.prins.err',
         gem='log/taxmaps/{seq}.{db}.gem.err',
         time='time/taxmaps/{seq}.{db}.log'
    params: db=taxmaps_db,
            tax=taxmaps_tax
    threads: 16
    resources: mem=60
    shell:
        taxmaps_rule

rule taxmaps_benchmark:
    input: fastq_both_input
    output: report='benchmark/reports/{seq}.taxmaps.{db}.tsv',
            summary='benchmark/extra_reports/{seq}.taxmaps.lca.summary.{db}.tsv',
            lca='benchmark/data/{seq}.taxmaps.lca.{db}.tsv'
    priority: 2
    log: prins='benchmark/log/taxmaps/{seq}.{db}.prins.err',
            gem='benchmark/log/taxmaps/{seq}.{db}.gem.err',
            time='benchmark/time/taxmaps/{seq}.{db}.log'
    params: db=taxmaps_db,
            tax=taxmaps_tax
    threads: ALL_CORES
    resources: mem=60
    benchmark: repeat('benchmark/{seq}/taxmaps.{db}.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('{DROPCACHE}')
        shell(taxmaps_rule, bench_record=bench_record)

TAXMAPS_TBL = config.get('TAXMAPS_TBL', 'taxMaps-taxtbl')
TAXMAPS_INDEX = config.get('TAXMAPS_INDEX', 'taxMaps-index')

rule taxmaps_refseqc_db:
    input: fna='db/refseqc/fasta/genomic.gi.fna'
    output: 'db/refseqc/taxmaps/refseqc.gem'
    params: dir='db/refseqc/taxmaps'
    log: log='log/db/taxmaps/refseqc.log',
         time='time/db/taxmaps/refseqc.log'
    benchmark: 'benchmark/db/taxmaps/refseqc.tsv'
    run:
        shell('{DROPCACHE}')
        shell('''\
        /usr/bin/time -v -o {log.time}  \
        env PATH={GEMPATH}:$PATH \
        {TAXMAPS_TBL} -n {TAXONOMY_DB}/names.dmp -t {TAXONOMY_DB}/nodes.dmp > {params.dir}/taxonomy.tbl 2>&1 | tee {log.log}
        /usr/bin/time -a -v -o {log.time}  \
        env PATH={GEMPATH}:$PATH \
        {TAXMAPS_INDEX} -i {input.fna} -c {TAXONOMY_DB}/gi_taxid_nucl.dmp -t taxonomy.tbl -p {params.dir}/refseqc | tee -a {log.log}
        ''', bench_record=bench_record)
