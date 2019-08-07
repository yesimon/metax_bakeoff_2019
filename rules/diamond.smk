import os
from os.path import join

DIAMOND = config.get('DIAMOND', 'diamond')

diamond_dbs = ['nr', 'refseqc']

DIAMOND_NR_ALL = expand('reports/{seq}.diamondkraken.nr.report', seq=samples_se + samples_pe)
rule diamond_nr_all:
    input: DIAMOND_NR_ALL

DIAMOND_REFSEQC_ALL = expand('reports/{seq}.diamondkraken.refseqc.report', seq=samples_se + samples_pe)
rule diamond_refseqc_all:
    input: DIAMOND_REFSEQC_ALL

DIAMOND_ALL = DIAMOND_NR_ALL + DIAMOND_REFSEQC_ALL
rule diamond_all:
    input: DIAMOND_NR_ALL + DIAMOND_REFSEQC_ALL

DIAMOND_FASTA_ALL = expand('reports/{seq}.diamondkraken.{db}.report', seq=samples_fasta, db=['nr', 'refseqc'])
rule diamond_fasta_all:
    input: DIAMOND_FASTA_ALL

DIAMOND_BENCHMARK_ALL = expand('benchmark/reports/{seq}.diamondkraken.refseqc.report', seq=benchmark_samples_se + benchmark_samples_pe)
rule diamond_benchmark_all:
    input: DIAMOND_BENCHMARK_ALL

def diamond_db_args(wildcards):
    args = {}
    if wildcards.db == 'nr':
        args['dmnd'] = config['DIAMOND_NR_DB']
    elif wildcards.db == 'refseqc':
        args['dmnd'] = config['DIAMOND_REFSEQC_DB']
    else:
        raise Exception
    db_dir = os.path.dirname(args['dmnd'])
    args['taxonmap'] = join(db_dir, 'prot.accession2taxid.gz')
    args['taxonnodes'] = join(db_dir, 'nodes.dmp')
    return args

def diamond_input(wildcards):
    if wildcards.seq in samples_pe:
        return expand('fastq/{seq}{pair}.fastq.bz2', pair=paired_suffix, seq=wildcards.seq)
    elif wildcards.seq in samples_se:
        return 'fastq/{}.fastq.bz2'.format(wildcards.seq)
    elif wildcards.seq in samples_fasta:
        return expand('fastq/{seq}.fasta', seq=wildcards.seq)
    else:
        raise Exception


# def diamond_input_args(wildcards, input):
#     if wildcards.seq in samples_pe:
#         return r'''\
#         paste <(paste - - - - < <(lbzcat {}) | awk 'BEGIN {{FS=OFS="\t"}} {{sub(/\/[1,2]$/, "", $1); print $1, $2, $3, $4 }}') \
#         <(paste - - - - < <(lbzcat {}) | awk 'BEGIN {{FS=OFS="\t"}} {{sub(/\/[1,2]$/, "", $1); print $1, $2, $3, $4 }}') \
#         | tr '\t' '\n' '''.format(input[0], input[1])
#     elif wildcards.seq in samples_se:
#         return 'lbzcat {}'.format(input)
#     elif wildcards.seq in samples_fasta:
#         return 'cat {}'.format(input)
#     else:
#         raise Exception

rule diamond:
    input: fastq_both_input
    # input: diamond_input
    output: tax='data/{seq}.diamond_tax.{db}.tsv.gz'
    params: dmnd=lambda wildcards: diamond_db_args(wildcards)['dmnd'],
            taxonmap=lambda wildcards: diamond_db_args(wildcards)['taxonmap'],
            taxonnodes=lambda wildcards: diamond_db_args(wildcards)['taxonnodes'],
            # input_args=diamond_input_args,
            blocksize=' -b ' + str(config['DIAMOND_BLOCKSIZE']) if 'DIAMOND_BLOCKSIZE' in config else '',
            index_chunks=' --index-chunks ' + str(config['DIAMOND_CHUNKS']) if 'DIAMOND_CHUNKS' in config else '',
            tmpdir=' --tmpdir ' + config['DIAMOND_TMPDIR'] if 'DIAMOND_TMPDIR' in config else '',
            tax='data/{seq}.diamond_tax.{db}.tsv'
    log: log='log/diamond/{seq}.{db}.log',
         time='time/diamond/{seq}.{db}.log'
    benchmark: 'benchmark/diamond/{seq}.{db}.log'
    threads: ALL_CORES
    shell:
        '''
        /usr/bin/time --verbose --append -o {log.time} \
        {DIAMOND} blastx -q {input} --verbose -p {threads} -d {params.dmnd}{params.tmpdir}{params.blocksize}{params.index_chunks} --taxonmap {params.taxonmap} --taxonnodes {params.taxonnodes} --outfmt 102 -o {params.tax} 2>&1 | tee {log.log}
        {PIGZ} -9 {params.tax}
        '''

rule diamond_benchmark:
    input: fastq_both_input
    output: tax='benchmark/data/{seq}.diamond_tax.{db}.tsv',
            report='benchmark/reports/{seq}.diamondkraken.{db}.report'
    params: dmnd=lambda wildcards: diamond_db_args(wildcards)['dmnd'],
            taxonmap=lambda wildcards: diamond_db_args(wildcards)['taxonmap'],
            taxonnodes=lambda wildcards: diamond_db_args(wildcards)['taxonnodes'],
            # input_args=diamond_input_args,
            blocksize=' -b ' + str(config['DIAMOND_BLOCKSIZE']) if 'DIAMOND_BLOCKSIZE' in config else '',
            index_chunks=' --index-chunks ' + str(config['DIAMOND_CHUNKS']) if 'DIAMOND_CHUNKS' in config else '',
            tmpdir=' --tmpdir ' + config['DIAMOND_TMPDIR'] if 'DIAMOND_TMPDIR' in config else ''
    log: log='benchmark/log/diamond/{seq}.{db}.log',
         time='benchmark/time/diamond/{seq}.{db}.log'
    benchmark: repeat('benchmark/{seq}/diamond.{db}.log', 2)
    threads: ALL_CORES
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell('''\
        {params.input_args} | \
        /usr/bin/time --verbose --append -o {log.time} \
        {DIAMOND} blastx --verbose --log -p {threads} -d {params.dmnd}{params.tmpdir}{params.blocksize}{params.index_chunks} --taxonmap {params.taxonmap} --taxonnodes {params.taxonnodes} --outfmt 102 -o {output.tax} 2>&1 | tee {log.log}

        ''', bench_record=bench_record)
        shell('''\
        metax taxid-report --output {output.report} --tax-dir {config[TAXONOMY_DB]} {output.tax}
        ''')

rule diamond_tax_lca:
    input: 'data/{seq}.diamond_tax.{db}.tsv.gz'
    output: 'reports/{seq}.diamondkraken.{db}.report'
    resources: io=50
    shell:
        '''
        metax taxid-report --output {output} --tax-dir {config[TAXONOMY_DB]} <({PIGZ} -dc {input})
        '''

rule diamond_refseqc_db:
    input: faa='db/refseqc/fasta/protein.uniq.faa.gz'
    output: 'db/refseqc/diamond/refseqc.dmnd'
    benchmark: 'benchmark/db/diamond/refseqc.tsv'
    run:
        shell('dropcache')
        shell('''\
        {DIAMOND} makedb --in {input.faa} -d db/refseqc/diamond/refseqc
        ''', bench_record=bench_record)
