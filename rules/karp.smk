import textwrap
KARP = config.get('KARP', 'karp')

def karp_input_args(wildcards, input):
    if wildcards.seq in samples_pe:
        return '-f {} -q {}'.format(input[0], input[1])
    elif wildcards.seq in samples_se:
        return '-f {}'.format(input)
    elif wildcards.seq in samples_fasta:
        return '-f {}'.format(input)
    else:
        raise Exception

rule karp_all:
    input: expand('reports/{sample}.karp.nofilter.collapse.freqs', sample=samples_all)
    # Karp seems to have problems with filter and nocollapse stringently removing too many reads
    # input: expand('reports/{sample}.karp.nofilter.freqs', sample=samples_all)

rule karp_benchmark_all:
    input: expand('benchmark/reports/{sample}.karp.nofilter.collapse.freqs', sample=benchmark_samples_all)

KARP_SHELL = '''\
{KARP} -c quantify --threads {threads} -r {config[KARP_SILVA_PREFIX]}.fasta -i {config[KARP_SILVA_PREFIX]}.index -t {config[KARP_SILVA_PREFIX]}.tax {params.input_args} --out {params.prefix} --like_thresh 19 --no_harp_filter --collapse
mv {params.prefix}.freqs {output}
cat {params.prefix}.log >> {log}
'''
rule karp_nofilter_collapse:
    input: fastq_input
    output: 'reports/{seq}.karp.nofilter.collapse.freqs'
    params: prefix='{seq}.karp',
            input_args=karp_input_args
    log: 'log/karp/{seq}.nofilter.collapse.log'
    threads: ALL_CORES
    resources: mem=100
    shadow: 'shallow'
    shell:
        KARP_SHELL


rule karp_nofilter_collapse_benchmark:
    input: fastq_input
    output: 'benchmark/reports/{seq}.karp.nofilter.collapse.freqs'
    params: prefix='{seq}.karp',
            input_args=karp_input_args
    log: 'benchmark/log/karp/{seq}.nofilter.collapse.log'
    threads: ALL_CORES
    resources: mem=100
    shadow: 'shallow'
    benchmark: repeat('benchmark/{seq}/karp.log', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(KARP_SHELL, bench_record=bench_record)

# rule karp:
#     input: karp_input
#     output: 'reports/{seq}.karp.freqs'
#     params: prefix='{seq}.karp',
#             input_args=karp_input_args
#     log: 'log/karp/{seq}.log'
#     threads: ALL_CORES
#     resources: mem=100
#     shadow: 'shallow'
#     shell:
#         '''
#         {KARP} -c quantify --threads {threads} -r {config[KARP_SILVA_PREFIX]}.fasta -i {config[KARP_SILVA_PREFIX]}.index -t {config[KARP_SILVA_PREFIX]}.tax {params.input_args} --out {params.prefix}
#         mv {params.prefix}.freqs {output}
#         cat {params.prefix}.log >> {log}
#         '''
