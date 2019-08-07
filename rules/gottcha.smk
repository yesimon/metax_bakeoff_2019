GOTTCHA = config.get('GOTTCHA')


GOTTCHA_SPECIES_ALL = expand('reports/{sample}.gottcha.tsv', sample=samples_all)
rule gottcha_species_all:
    input: GOTTCHA_SPECIES_ALL

GOTTCHA_GENUS_ALL = expand('reports/{sample}.gottcha_genus.tsv', sample=samples_all)
rule gottcha_genus_all:
    input: GOTTCHA_GENUS_ALL

GOTTCHA_ALL = expand('reports/{sample}.gottcha.tsv', sample=samples_all)
rule gottcha_all:
    input: GOTTCHA_ALL

GOTTCHA_BENCHMARK_ALL = expand('benchmark/reports/{sample}.gottcha.tsv', sample=benchmark_samples_all)
rule gottcha_benchmark_all:
    input: GOTTCHA_BENCHMARK_ALL

def gottcha_bacteria_db_args(wildcards, input, output):
    if '_genus' not in output[0]:
        return config['GOTTCHA_BACTERIA_DB']
    else:
        return config['GOTTCHA_BACTERIA_DB'].replace('species', 'genus')

def gottcha_virus_db_args(wildcards, input, output):
    if '_genus' not in output[0]:
        return config['GOTTCHA_VIRUS_DB']
    else:
        return config['GOTTCHA_VIRUS_DB'].replace('species', 'genus')

GOTTCHA_SHELL = '''\
{GOTTCHA} --database {params.bacteria_db} --input {input} --mode all --dumpSam --prefix {wildcards.seq} --threads {threads}
mv {wildcards.seq}.gottcha_full.tsv {wildcards.seq}.bacteria.tsv || touch {wildcards.seq}.bacteria.tsv
touch -a {wildcards.seq}.gottcha.sam
samtools view -h -b --threads {threads} -o {output.bacteria_sam} {wildcards.seq}.gottcha.sam
{GOTTCHA} --database {params.virus_db} --input {input} --mode all --dumpSam --prefix {wildcards.seq} --threads {threads}
mv {wildcards.seq}.gottcha_full.tsv {wildcards.seq}.virus.tsv || touch {wildcards.seq}.virus.tsv
touch -a {wildcards.seq}.gottcha.sam
samtools view -h -b --threads {threads} -o {output.virus_sam} {wildcards.seq}.gottcha.sam

sed -i '1{{/LEVEL/d}}' {wildcards.seq}.virus.tsv

cat {wildcards.seq}.bacteria.tsv {wildcards.seq}.virus.tsv > {output.report} || true
'''
rule gottcha:
    input: fastq_both_input
    output: report='reports/{seq}.gottcha{rank,(_genus)?}.tsv',
            bacteria_sam='data/{seq}.gottcha{rank,(_genus)?}.bacteria.bam',
            virus_sam='data/{seq}.gottcha{rank,(_genus)?}.virus.bam'
    params: qualfix='{seq}.qualfix.fastq',
            bacteria_db=gottcha_bacteria_db_args,
            virus_db=gottcha_virus_db_args
    priority: 1
    shadow: 'shallow'
    threads: ALL_CORES
    resources: mem=16
    shell:
        GOTTCHA_SHELL

rule gottcha_benchmark:
    input: fastq_both_input
    output: report='benchmark/reports/{seq}.gottcha.tsv',
            bacteria_sam='benchmark/data/{seq}.gottcha.bacteria.bam',
            virus_sam='benchmark/data/{seq}.gottcha.virus.bam'
    params: qualfix='{seq}.qualfix.fastq',
            bacteria_db=gottcha_bacteria_db_args,
            virus_db=gottcha_virus_db_args
    priority: 1
    shadow: 'shallow'
    threads: ALL_CORES
    resources: mem=16
    benchmark: repeat('benchmark/{seq}/gottcha.tsv', 2)
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(GOTTCHA_SHELL, bench_record=bench_record)
        shell('truncate -s 0 {output}')
