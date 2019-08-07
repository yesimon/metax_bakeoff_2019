PHYLOSIFT = config.get('PHYLOSIFT')


def phylosift_input(wildcards):
    if wildcards.seq in samples:
        return expand('fastq/{}.{{pair}}.fastq'.format(wildcards.seq), pair=[1, 2])
    else:
        return 'fastq/{}.fastq'.format(wildcards.seq)

def phylosift_paired(wildcards, input):
    if len(input) > 1:
        return ' --paired'
    return ''

rule phylosift_all:
    input: expand('reports/{sample}.phylosift.tsv', sample=samples_se + samples_pe)

rule phylosift:
    input: phylosift_input
    output: data='data/{seq}.phylosift.tsv.gz',
            data_lineage='data/{seq}.phylosift_lineage.tsv.gz',
            report='reports/{seq}.phylosift.tsv'
    params: output_dir='tmp/phylosift/{seq}',
            paired=phylosift_paired
    shadow: 'shallow'
    log: log='log/phylosift/{seq}.log',
         time='time/phylosift/{seq}.log'
    threads: ALL_CORES
    shell:
       '''
       /usr/bin/time -v -o {log.time} \
         {PHYLOSIFT} all --output={params.output_dir}{params.paired} {input} | tee {log.log}
       {PIGZ} -c {params.output_dir}/sequence_taxa.txt > {output.data}
       {PIGZ} -c {params.output_dir}/sequence_taxa_summary.txt > {output.data_lineage}
       mv {params.output_dir}/taxasummary.txt {output.report}
       '''
