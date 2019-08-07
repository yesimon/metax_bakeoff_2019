

def motus_input(wildcards, input):
    if wildcards.seq in samples_fasta:
        tpl = ' -s {}'.format(input)
    else:
        if wildcards.seq in samples_pe:
            tpl = ' -f {} -r {}'.format(input[0], input[1])
        elif wildcards.seq in samples_se:
            tpl = ' -s {}'.format(input)
        else:
            raise Exception
    return tpl

MOTUS_ALL = expand('reports/{seq}.motus.txt', seq=samples_all)
rule motus2_all:
    input: MOTUS_ALL


MOTUS_BENCHMARK_ALL = expand('benchmark/reports/{seq}.motus.txt', seq=benchmark_samples_all)
rule motus2_benchmark_all:
    input: MOTUS_BENCHMARK_ALL

MOTUS_SHELL = textwrap.dedent('''\
/usr/bin/time -v -o {log.time} \
motus profile -puc {params.input} -o {output.report} -t {threads} 2>&1 | tee {log.log}
''')

rule motus2:
    input: fastq_input
    output: report='reports/{seq}.motus.txt'
    params: input=motus_input
    log: log='log/motus2/{seq}.log',
         time='time/motus2/{seq}.log'
    threads: ALL_CORES
    resources: io=50
    shell:
        MOTUS_SHELL

rule motus2_benchmark:
    input: fastq_input
    output: report='benchmark/reports/{seq}.motus.txt'
    params: input=motus_input
    log: log='benchmark/log/motus2/{seq}.log',
         time='benchmark/time/motus2/{seq}.log'
    threads: ALL_CORES
    resources: io=50
    benchmark: repeat('benchmark/{seq}/motus2.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(MOTUS_SHELL, bench_record=bench_record)
        shell('truncate -s 0 {output}')
