import os
from os.path import join

clark_base_all = expand('reports/{sample}.{{exe}}{{rank}}.{{db}}.csv', sample=samples_all)
CLARK_ALL = expand(clark_base_all, exe=['clark', 'clark_s'], rank=['', '_genus'], db=['default_bv', 'refseqc'])
rule clark_all:
    input: CLARK_ALL

rule clark_default_species_all:
    input: expand(clark_base_all, exe='clark', rank='', db='default_bv')

rule clark_default_genus_all:
    input: expand(clark_base_all, exe='clark', rank='_genus', db='default_bv')

rule clark_s_default_species_all:
    input: expand(clark_base_all, exe='clark_s', rank='', db='default_bv')

rule clark_s_default_genus_all:
    input: expand(clark_base_all, exe='clark_s', rank='_genus', db='default_bv')

rule clark_refseqc_species_all:
    input: expand(clark_base_all, exe='clark', rank='', db='refseqc')

rule clark_refseqc_genus_all:
    input: expand(clark_base_all, exe='clark', rank='_genus', db='refseqc')

rule clark_s_refseqc_species_all:
    input: expand(clark_base_all, exe='clark_s', rank='', db='refseqc')

rule clark_s_refseqc_genus_all:
    input: expand(clark_base_all, exe='clark_s', rank='_genus', db='refseqc')

CLARK_BENCHMARK_ALL = expand('benchmark/reports/{sample}.clark.refseqc.csv', sample=benchmark_samples_all)
rule clark_benchmark_all:
    input: CLARK_BENCHMARK_ALL

CLARK_S_BENCHMARK_ALL = expand('benchmark/reports/{sample}.clark_s.refseqc.csv', sample=benchmark_samples_all)
rule clark_s_benchmark_all:
    input: CLARK_S_BENCHMARK_ALL

def clark_input_args(wildcards, input):
    if wildcards.seq in samples_pe:
        return '-P {} {}'.format(input[0], input[1])
    elif wildcards.seq in samples_se:
        return '-O {}'.format(input)
    elif wildcards.seq in samples_fasta:
        return '-O {}'.format(input)
    else:
        raise Exception

def clark_db_args(wildcards, input, output):
    if output[0].endswith('default_bv.csv.gz'):
        db = config.get('CLARK_BV_DB')
    elif output[0].endswith('refseqc.csv.gz'):
        db = config.get('CLARK_REFSEQC_DB')
    else:
        raise Exception

    if '_genus' in output[0]:
        db = db + '_1'
    else:
        db = db + '_0'

    targets = join(db, 'targets.txt')
    return '-D {} -T {}'.format(db, targets)

def clark_exe(wildcards, input, output):
    if wildcards.exe == 'clark':
        return config.get('CLARK', 'CLARK')
    elif wildcards.exe == 'clark_s':
        return config.get('CLARK_S', 'CLARK-S')

def clark_mem(wildcards):
    if wildcards.db == 'refseqc':
        return 100
    elif wildcards.db == 'default_bv':
        return 190
    else:
        raise Exception

CLARK_SHELL = '''\
/usr/bin/time -v -o {log.time} \
{params.exe} -n {threads} {params.db} {params.input_args} -R {params.output_prefix} 2>&1 | tee {log.log}
{PIGZ} -9 {params.data}
'''

rule clark:
    input: fastx_input
    output: data='data/{seq}.{exe}.{db}.csv.gz'
    params: exe=clark_exe,
            db=clark_db_args,
            input_args=clark_input_args,
            data='data/{seq}.{exe}.{db}.csv',
            output_prefix='data/{seq}.{exe}.{db}'
    wildcard_constraints:
        exe="(clark)|(clark_s)"
    priority: 3
    log: log='log/{exe}/{seq}.{db}.log',
         time='time/{exe}/{seq}.{db}.log'
    threads: ALL_CORES
    resources: mem=160
    shell:
        CLARK_SHELL

rule clark_genus:
    input: fastx_input
    output: data='data/{seq}.{exe}_genus.{db}.csv.gz'
    params: exe=clark_exe,
            db=clark_db_args,
            input_args=clark_input_args,
            data='data/{seq}.{exe}_genus.{db}.csv',
            output_prefix='data/{seq}.{exe}_genus.{db}'
    wildcard_constraints:
        exe="(clark)|(clark_s)"
    priority: 3
    log: log='log/{exe}_genus/{seq}.{db}.log',
         time='time/{exe}_genus/{seq}.{db}.log'
    resources: mem=clark_mem
    threads: ALL_CORES
    shell:
        CLARK_SHELL

def clark_abundance_db_args(wildcards, input, output):
    if output[0].endswith('default_bv.csv'):
        db = config.get('CLARK_BV_DB')
    elif output[0].endswith('refseqc.csv'):
        db = config.get('CLARK_REFSEQC_DB')
    else:
        raise Exception

    if '_genus' in output[0]:
        db = db + '_1'
    else:
        db = db + '_0'

    return '-D {}'.format(os.path.dirname(db))

def abundance_clark_exe(wildcards, input, output):
    dirname = clark_exe(wildcards, input, output)
    return os.path.join(os.path.dirname(dirname), 'getAbundance')

rule clark_abundance:
    input: 'data/{seq}.{exe}{rank,(_genus)?}.{db}.csv.gz'
    output: 'reports/{seq}.{exe}{rank,(_genus)?}.{db}.csv'
    params: exe=abundance_clark_exe,
            db=clark_abundance_db_args
    wildcard_constraints:
        exe="(clark)|(clark_s)"
    shell:
        '''
        {params.exe} -F <(zcat {input}) {params.db} > {output}
        '''

CLARK_BENCHMARK_ABUNDANCE_SHELL = '{params.abundance_exe} -F {output.data} {params.abundance_db} > {output.report}'
rule clark_benchmark:
    input: fastx_input
    output: data='benchmark/data/{seq}.{exe}.{db}.csv',
            report='benchmark/reports/{seq}.{exe}.{db}.csv'
    params: exe=clark_exe,
            abundance_exe=abundance_clark_exe,
            db=clark_db_args,
            abundance_db=clark_abundance_db_args,
            input_args=clark_input_args,
            output_prefix='benchmark/data/{seq}.{exe}.{db}'
    wildcard_constraints:
        exe="(clark)|(clark_s)"
    priority: 3
    log: log='benchmark/log/{exe}/{seq}.{db}.log',
         time='benchmark/time/{exe}/{seq}.{db}.log'
    threads: ALL_CORES
    resources: mem=clark_mem
    benchmark: repeat('benchmark/{seq}/{exe}.{db}.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('{DROPCACHE}')
        shell(CLARK_SHELL + CLARK_BENCHMARK_ABUNDANCE_SHELL, bench_record=bench_record)
        shell('truncate -s 0 {output}')

rule clark_refseqc_db:
    input: fna='db/refseqc/fasta/genomic.gi.fna'
    output: 'db/refseqc/clark/refseqc/custom_0/db_central_k31_t10675_s1610612741_m0.tsk.ky'
    params: dir='db/refseqc/clark'
    log: log='log/db/clark/refseqc.log',
         time='time/db/clark/refseqc.log'
    benchmark: 'benchmark/db/clark/refseqc.tsv'
    run:
        shell('{DROPCACHE}')
        shell('''\
        /usr/bin/time -v -o {log.time}  \
        {config[CLARK_OPT]}/set_targets.sh /mnt/metax/workflow/db/refseqc/clark custom
        {config[CLARK]} -D {params.dir}/custom_0 -T {params.dir}/targets.txt -O /mnt/metax/workflow/tmp/empty.fastq -R empty.results
        ''', bench_record=bench_record)

rule clark_s_refseqc_db:
    input: 'db/refseqc/clark/custom_0/db_central_k31_t10673_s1610612741_m0.tsk.ky'
    output: 'db/refseqc/clark/custom_0/T58570/db_central_k31_t10673_s1610612741_m0_w22.tsk.ky'
    params: dir='db/refseqc/clark'
    log: log='log/db/clark_s/refseqc.log',
         time='time/db/clark_s/refseqc.log'
    benchmark: 'benchmark/db/clark_s/refseqc.tsv'
    run:
        shell('{DROPCACHE}')
        shell('''\
        echo "$(readlink -f {params.dir}/custom_0)" > {config[CLARK_OPT]}/.dbAddress
        # Must be in the CLARK install directory
        cd {config[CLARK_OPT]}
        {config[CLARK_OPT]}/buildSpacedDB.sh
        ''', bench_record=bench_record)
