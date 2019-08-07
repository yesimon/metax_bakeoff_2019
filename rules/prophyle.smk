PROPHYLE = config.get('PROPHYLE', 'prophyle')

PROPHYLE_REFSEQC_ALL = expand('reports/{sample}.prophyle.refseqc.txt', sample=samples_all)
PROPHYLE_BVP_ALL = expand('reports/{sample}.prophyle.bvp.txt', sample=samples_all)
PROPHYLE_ALL = PROPHYLE_REFSEQC_ALL + PROPHYLE_BVP_ALL
rule prophyle_all:
    input: PROPHYLE_ALL

rule prophyle_default_all:
    input: PROPHYLE_BVP_ALL

rule prophyle_refseqc_all:
    input: PROPHYLE_REFSEQC_ALL

PROPHYLE_FASTA_ALL = expand('reports/{sample}.prophyle.{db}.txt', sample=samples_fasta, db=['refseqc', 'bvp'])
rule prophyle_fasta_all:
    input: PROPHYLE_FASTA_ALL

PROPHYLE_BENCHMARK_ALL = expand('benchmark/data/{sample}.prophyle.refseqc.sam', sample=samples_all)
rule prophyle_benchmark_all:
    input: PROPHYLE_BENCHMARK_ALL

def prophyle_db(wildcards):
    if wildcards.db == 'refseqc':
        return config['PROPHYLE_REFSEQC_DB']
    elif wildcards.db == 'bvp':
        return config['PROPHYLE_BVP_DB']
    else:
        raise Exception

def prophyle_paired(wildcards):
    if wildcards.seq in samples_pe:
        return ' --paired'
    else:
        return ''

PROPHYLE_SHELL = '''\
/usr/bin/time -v -o {log.time} \
{PROPHYLE} classify {params.db} {input} > {output.data} 2> {log.log}
'''
rule prophyle:
    input: fastx_input
    output: data='data/{seq}.prophyle.{db}.sam'
    params: db=prophyle_db
    log: log='log/prophyle/{seq}.{db}.log',
         time='time/prophyle/{seq}.{db}.log'
    benchmark: 'benchmark/prophyle/{seq}.{db}.log'
    threads: 1
    resources: mem=30, io=10
    shell:
        PROPHYLE_SHELL

rule prophyle_benchmark:
    input: fastx_input
    output: data='benchmark/data/{seq}.prophyle.{db}.sam'
    params: db=prophyle_db
    log: log='benchmark/log/prophyle/{seq}.{db}.log',
            time='benchmark/time/prophyle/{seq}.{db}.log'
    benchmark: repeat('benchmark/{seq}/prophyle.{db}.log', 2)
    threads: ALL_CORES
    resources: mem=30
    run:
        if benchmark_i == 0:
            shell('dropcache')
        shell(PROPHYLE_SHELL)

rule prophyle_lca:
    input: data='data/{seq}.prophyle.{db}.sam',
           total_reads='info/{seq}.total_reads.txt'
    output: report='reports/{seq}.prophyle.{db}.txt'
    params: db=TAXONOMY_DB,
            paired=prophyle_paired
    shell:
        '''
        metax prophyle-report{params.paired} --total-reads $(cat {input.total_reads}) --output {output.report} --tax-dir {params.db} {input.data}
        '''

rule prophyle_refseqc_db:
    input: assembly_summary=join(NCBI_GENOMES, 'assembly_summary.txt')
    params: tax_db='db/taxonomy/20180425'
    log: log='log/db/prophyle/refseqc.log',
         time='time/db/prophyle/refseqc.log'
    benchmark: 'benchmark/db/prophyle/refseqc.tsv'
    threads: ALL_CORES
    resources: mem=60
    shadow: 'shallow'
    run:
        shell('''
        mkdir -p db/refseqc/build_prophyle

        python -c "from ete3 import ncbi_taxonomy; ncbi_taxonomy.NCBITaxa(taxdump_file='{params.tax_db}/taxdump.tar.gz')"
        dropcache
        ''')

        shell('''\
        awk -F "\t" -v OFS="\t" '$12=="Complete Genome" && $11=="latest" {{print $1, $6, $20}}' {input.assembly_summary} > ftpselection.tsv
        cut -f 3 ftpselection.tsv | awk 'BEGIN{{FS=OFS="/";filesuffix="genomic.fna.gz"}} {{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}}' > ftpfilepaths.tsv
        cut -f 1,2 ftpselection.tsv | sed 's/\.[0-9]*//' > acc2taxid.tsv

        prophyle_ncbi_tree.py archaea {NCBI_GENOMES}/refseq archaea.nw acc2taxid.tsv -l archaea.log
        prophyle_ncbi_tree.py bacteria {NCBI_GENOMES}/refseq bacteria.nw acc2taxid.tsv -l bacteria.log
        prophyle_ncbi_tree.py viral {NCBI_GENOMES}/refseq viral.nw acc2taxid.tsv -l viral.log

        ln -s {NCBI_GENOMES}/refseq/archaea .
        ln -s {NCBI_GENOMES}/refseq/bacteria .
        ln -s {NCBI_GENOMES}/refseq/viral .
        /usr/bin/time -v -o {log.time} \
          prophyle index -k 31 archaea.nw bacteria.nw viral.nw refseqc | tee {log.log}
        ''' , bench_record=bench_record)
        shell('mv refseqc ../prophyle')
