AC_DIAMOND = config.get('AC_DIAMOND', 'ac-diamond')

AC_DIAMOND_NR_ALL = expand('data/{seq}.ac_diamond.{db}.m8.zst', seq=samples_se + samples_pe, db=['nr'])
AC_DIAMOND_REFSEQC_ALL = expand('data/{seq}.ac_diamond.{db}.m8.zst', seq=samples_se + samples_pe, db=['refseqc'])
# AC_DIAMOND_ALL = expand('reports/{seq}.diamondkraken.{db}.report', seq=samples_se + samples_pe, db=['refseqc'])
# AC_DIAMOND_ALL = expand('data/{seq}.ac_diamond.{db}.m8.zst', seq=samples_se + samples_pe, db=['refseqc'])
rule ac_diamond_all:
    input: AC_DIAMOND_NR_ALL + AC_DIAMOND_REFSEQC_ALL

rule ac_diamond_report_all:
    input: expand('reports/{seq}.ac_diamond.{db}.txt', seq=samples_se + samples_pe, db=['refseqc'])

rule ac_diamond_refseqc_all:
    input: AC_DIAMOND_REFSEQC_ALL

rule ac_diamond_nr_all:
    input: AC_DIAMOND_NR_ALL

AC_DIAMOND_SENSITIVE_ALL = expand('data/{seq}.ac_diamond_sensitive.{db}.daa', seq=samples_se + samples_pe, db=['refseqc'])
rule ac_diamond_sensitive_all:
    input: AC_DIAMOND_SENSITIVE_ALL

AC_DIAMOND_FASTA_ALL = expand('reports/{seq}.diamondkraken.{db}.report', seq=samples_fasta, db=['refseqc'])
rule ac_diamond_fasta_all:
    input: AC_DIAMOND_FASTA_ALL

AC_DIAMOND_BENCHMARK_ALL = expand('benchmark/data/{seq}.ac_diamond.refseqc.daa', seq=benchmark_samples_se + benchmark_samples_pe)
rule ac_diamond_benchmark_all:
    input: AC_DIAMOND_BENCHMARK_ALL

def ac_diamond_db_args(wildcards):
    if wildcards.sens:
        if wildcards.db == 'nr':
            return config['AC_DIAMOND_NR_SENSITIVE_DB']
        elif wildcards.db == 'refseqc':
            return config['AC_DIAMOND_REFSEQC_SENSITIVE_DB']
    else:
        if wildcards.db == 'nr':
            return config['AC_DIAMOND_NR_DB']
        elif wildcards.db == 'refseqc':
            return config['AC_DIAMOND_REFSEQC_DB']


def ac_diamond_input_args(wildcards):
    if wildcards.sens:
        return ' --sensitive'
    else:
        return ''

AC_DIAMOND_SHELL = textwrap.dedent('''\
/usr/bin/time --verbose --append -o {log.time} \
{AC_DIAMOND} align --verbose -q {input} -p {threads} -a {output.daa} -d {params.dmnd}{params.tmpdir}{params.input_args} 2>&1 | tee {log.log}
{AC_DIAMOND} view --verbose -p {threads} -a {output.daa} -o {params.m8} 2>&1 | tee {log.log}

zstd -19 -T{threads} {params.m8}
rm {params.m8}
''')


rule ac_diamond:
    input: fastq_both_input
    output: daa='data/{seq}.ac_diamond{sens,(_sensitive)?}.{db}.daa',
            m8='data/{seq}.ac_diamond{sens,(_sensitive)?}.{db}.m8.zst'
    params: dmnd=ac_diamond_db_args,
            input_args=ac_diamond_input_args,
            m8='data/{seq}.ac_diamond{sens}.{db}.m8',
            tmpdir=' --tmpdir ' + config['DIAMOND_TMPDIR'] if 'DIAMOND_TMPDIR' in config else ''
    log: log='log/ac_diamond/{seq}{sens}.{db}.log',
         time='time/ac_diamond/{seq}{sens}.{db}.log'
    benchmark: 'benchmark/ac_diamond/{seq}{sens}.{db}.log'
    threads: ALL_CORES
    shell:
        AC_DIAMOND_SHELL

rule ac_diamond_benchmark:
    input: fastq_both_input
    output: daa='benchmark/data/{seq}.ac_diamond{sens,(_sensitive)?}.{db}.daa',
            m8='benchmark/data/{seq}.ac_diamond{sens,(_sensitive)?}.{db}.m8.zst'
    params: dmnd=ac_diamond_db_args,
            input_args=ac_diamond_input_args,
            m8='benchmark/data/{seq}.ac_diamond{sens}.{db}.m8',
            tmpdir=' --tmpdir ' + config['DIAMOND_TMPDIR'] if 'DIAMOND_TMPDIR' in config else ''
    log: log='benchmark/log/ac_diamond/{seq}{sens}.{db}.log',
         time='benchmark/time/ac_diamond/{seq}{sens}.{db}.log'
    benchmark: repeat('benchmark/{seq}/ac_diamond{sens}.{db}.tsv', 2)
    threads: ALL_CORES
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(AC_DIAMOND_SHELL, bench_record=bench_record)
        shell('truncate -s 0 {output}')

rule ac_diamond_report:
    input: data='data/{seq}.ac_diamond.{db}.m8.zst',
           total_reads='info/{seq}.total_reads.txt'
    output: report='reports/{seq}.ac_diamond.{db}.txt'
    params: db=config.get('TAXONOMY_DB')
    shell:
        '''
        export TMPDIR={TMPDIR}
        metax blast-report --total-reads $(cat {input.total_reads}) --blast-report {output.report} --tax-dir {params.db} <(zstd -dc {input.data})
        '''
