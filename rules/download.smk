from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
GS = GSRemoteProvider()

db_dirs = ['blast/default',
           'blast/refseqc',
           'bracken/default',
           'bracken/refseqc',
           'clark/refseqc',
           'clark/default',
           'centrifuge/refseqc',
           'centrifuge/p_compressed+h+v',
           'diamond/nr',
           'diamond/refseqc',
           'kaiju/nr',
           'kaiju/refseqc',
           'karp/silva',
           'kraken/default',
           'kraken/refseqc',
           'kslam/default',
           'kslam/refseqc',
           'metaothello/default',
           'metaphlan2/default',
           'mmseqs2/nr',
           'mmseqs2/refseqc',
           'pathseq/default',
           'prophyle/default',
           'prophyle/refseqc',
           'taxmaps/default',
           'taxmaps/refseqc']

rule download_db_all:
    input: expand('db/{db}', db=db_dirs) + ['db/taxonomy']


rule download_taxonomy_db:
    input: GS.remote('metax-bakeoff-2019/db/taxdump.tar.gz', stay_on_remote=True)
    output: directory('db/taxonomy')
    shell:
        '''
        mkdir {output}
        gsutil cp {input} - | pigz -d | tar x -C {output}
        '''


rule download_db:
    input: GS.remote('metax-bakeoff-2019/db/{method}.{db}.tar.zst', stay_on_remote=True)
    output: directory('db/{method}/{db}')
    shell:
        '''
        mkdir {output}
        gsutil cp {input} - | zstd -d | tar x -C {output}
        '''


fastq_remote = ['{}_1'.format(x) for x in samples_pe] + ['{}_2'.format(x) for x in samples_pe] + samples_se



rule download_fastq_all:
    input: expand('fastq/{seq}.bz2', seq=fastq_remote)


rule download_fastq:
    input: GS.remote('metax-bakeoff-2019/fastq/{seq}.fastq.bz2', stay_on_remote=True)
    output: 'fastq/{seq}.fastq.bz2'
    shell:
        '''
        gsutil cp {input} {output}
        '''

rule download_src:
    input: GS.remote('metax-bakeoff-2019/src/src.tar.zst', stay_on_remote=True)
    output: directory('src')
    shell:
        '''
        gsutil cp {input} - | zstd -d | tar x
        '''
