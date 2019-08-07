BLASTN = config.get('BLASTN', 'blastn')
METAFLOW_DIR = config.get('METAFLOW_DIR')
METAFLOW_DB = config.get('METAFLOW_DB')

rule metaflow_all:
    input: expand('reports/{sample}.metaflow.tsv', sample=samples)

rule metaflow:
    input: expand('fastq/{{seq}}{pair}.fastq.gz', pair=paired_suffix)
    output: report='reports/{seq}.metaflow.tsv',
            blast='data/{seq}.metaflow_blast.m8'
    log: 'log/metaflow/{seq}.log'
    shadow: 'shallow'
    shell:
        '''
        {BLASTN} -query <(zcat {input} | seqtk seq -A - ) -out {wildcards.seq}.blast -outfmt 6 -db {METAFLOW_DB}/NCBI_DB/BLAST_DB.fasta
        python {METAFLOW_DIR}/BLAST_TO_LGF.py {wildcards.seq}.blast {METAFLOW_DB}/NCBI_DB/NCBI_Ref_Genome.txt 250 1 0.0001
        {METAFLOW_DIR}/metaflow -m {wildcards.seq}.blast.lgf -g {METAFLOW_DB}/NCBI_DB/NCBI_Ref_Genome.txt -c {METAFLOW_DIR}/metaflow.config
        mv {wildcards.seq}.blast {output.blast}
        mv {wildcards.seq}.blast.lgf.abundance.csv {output.report}
        '''


rule metaflow_single:
    input: 'fastq/{seq}.fastq.gz'
    output: report='reports/{seq}.metaflow.tsv',
            blast='data/{seq}.metaflow_blast.m8'
    log: 'log/metaflow/{seq}.log'
    shadow: 'shallow'
    shell:
        '''
        {BLASTN} -query <(zcat {input} | seqtk seq -A - ) -out {wildcards.seq}.blast -outfmt 6 -db {METAFLOW_DB}/NCBI_DB/BLAST_DB.fasta
        python {METAFLOW_DIR}/BLAST_TO_LGF.py {wildcards.seq}.blast {METAFLOW_DB}/NCBI_DB/NCBI_Ref_Genome.txt 250 1 0.0001
        {METAFLOW_DIR}/metaflow -m {wildcards.seq}.blast.lgf -g {METAFLOW_DB}/NCBI_DB/NCBI_Ref_Genome.txt -c {METAFLOW_DIR}/metaflow.config
        mv {wildcards.seq}.blast {output.blast}
        mv {wildcards.seq}.blast.lgf.abundance.csv {output.report}
        '''
