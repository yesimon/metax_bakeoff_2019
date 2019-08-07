CENTRIFUGE = config.get('CENTRIFUGE', 'centrifuge')


CENTRIFUGE_ALL = expand('reports/{sample}.centkraken.{db}.tsv', sample=samples_pe + samples_se, db=['refseqc', 'p_compressed+h+v'])
rule centrifuge_all:
    input: CENTRIFUGE_ALL

CENTRIFUGE_BENCHMARK_ALL = expand('benchmark/reports/{sample}.centkraken.refseqc.tsv', sample=benchmark_samples_pe + benchmark_samples_se)
rule centrifuge_benchmark_all:
    input: CENTRIFUGE_BENCHMARK_ALL

def centrifuge_input(wildcards, input):
    if wildcards.seq in samples_fasta:
        tpl = ' -f -U {}'.format(input)
    else:
        if wildcards.seq in samples_pe:
            tpl = ' -1 {} -2 {}'.format(input[0], input[1])
        elif wildcards.seq in samples_se:
            tpl = ' -U {}'.format(input)
        else:
            raise Exception
    return tpl

def centrifuge_db(wildcards):
    if wildcards.db == 'refseqc':
        return config.get('CENTRIFUGE_REFSEQC_DB')
    elif wildcards.db == 'p_compressed+h+v':
        return config.get('CENTRIFUGE_COMPRESSED_DB')
    else:
        raise Exception

CENTRIFUGE_SHELL = '''\
/usr/bin/time -v -o {log.time} \
{CENTRIFUGE} -x {params.db}{params.inf} --report-file {output.report} --threads {threads} | pigz -9 -f > {output.data} 2> {log.log}
'''
rule centrifuge:
    input: fastx_bz2_input
    output: report='reports/extra/{seq}.centrifuge.{db}.tsv',
            data='data/{seq}.centrifuge.{db}.tsv.gz',
    params: db=centrifuge_db,
            inf=centrifuge_input
    log: log='log/centrifuge/{db}/{seq}.log',
         time='time/centrifuge/{db}/{seq}.log'
    threads: ALL_CORES
    resources: mem=100
    shell:
        CENTRIFUGE_SHELL

rule centrifuge_kraken_report:
    input: 'data/{seq}.centrifuge.{db}.tsv.gz'
    output: report='reports/{seq}.centkraken.{db}.tsv'
    params: db=centrifuge_db
    shell:
        '''
        /mnt/metax/src/centrifuge/centrifuge-kreport -x {params.db} <(pigz -dc {input}) > {output}
        '''

CENTRIFUGE_REPORT_SHELL = '''\
/mnt/metax/src/centrifuge/centrifuge-kreport -x {params.db} <(pigz -dc {output.data}) > {output.kreport}
'''
rule centrifuge_bench:
    input: fastx_bz2_input
    output: report='benchmark/reports/{seq}.centrifuge.{db}.tsv',
            data='benchmark/data/{seq}.centrifuge.{db}.tsv.gz',
            kreport='benchmark/reports/{seq}.centkraken.{db}.tsv'
    params: db=centrifuge_db,
            inf=centrifuge_input
    log: log='benchmark/log/centrifuge/{db}/{seq}.log',
         time='benchmark/time/centrifuge/{db}/{seq}.log'
    benchmark: repeat('benchmark/{seq}/centrifuge.{db}.tsv', 2)
    threads: ALL_CORES
    resources: mem=100
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(CENTRIFUGE_SHELL + CENTRIFUGE_REPORT_SHELL, bench_record=bench_record)

rule centrifuge_refseqc_db:
    input: fna='db/refseqc/fasta/genomic.fna'
    output: 'db/refseqc/centrifuge/refseqc.1.cf'
    params: tax_db='db/taxonomy/20180425',
            dir='db/refseqc/centrifuge'
    log: log='log/db/centrifuge/refseqc.log',
         time='time/db/centrifuge/refseqc.log'
    benchmark: 'benchmark/db/centrifuge/refseqc.tsv'
    threads: ALL_CORES
    resources: mem=100
    run:
        # shell('''\
        # metax create-centrifuge-map {NCBI_GENOMES}

        # cat all-*.map > all.map
        # ''')
        shell('''\
        dropcache
        ''')
        shell('''\
        export TMPDIR={TMPDIR}
        TMP_DIR=$(mktemp -d)

	    /usr/bin/time -v -o {log.time} \
        centrifuge-build -p {threads} \
		  --conversion-table all.map \
		  --taxonomy-tree {params.tax_db}/nodes.dmp --name-table {params.tax_db}/names.dmp \
		  {input.fna} $TMP_DIR/refseqc 2>&1 | tee {log.log}

	    mv $TMP_DIR/refseqc.*.cf db/refseqc/centrifuge/ && rmdir $TMP_DIR
        ''', bench_record=bench_record)
