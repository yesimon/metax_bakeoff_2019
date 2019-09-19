
rule total_reads_all:
    input: expand('info/{seq}.total_reads.txt', seq=samples_all)

# Counts the number of reads (paired reads count as 2 reads)
rule total_reads:
    input: fastx_bz2_input
    output: 'info/{seq}.total_reads.txt'
    run:
        shell('echo {input}')
        if wildcards.seq in samples_fastq:
            shell('''\
            TOTAL_READS=$(expr $(wc -l <(lbzcat {input}) | cut -f1 -d' ') / 4)
            echo $TOTAL_READS > {output}
            ''')
        elif wildcards.seq in samples_fasta:
            shell('''\
            TOTAL_READS=$(fasta_formatter -i <(lbzcat {input}) -t | wc -l)
            echo $TOTAL_READS > {output}
            ''')


rule compile_total_reads:
    input: expand('info/{seq}.total_reads.txt', seq=samples_all)
    output: 'plotting/read_counts.tsv'
    shell:
        '''
        for i in {input}; do
            NAME="$(echo $(basename $i) | sed 's/.total_reads.txt//' )"
            echo -e "$NAME\t$(cat $i)" >> {output}
        done
        '''

rule compile_benchmark_all:
    input: expand('benchmark/{seq}.tsv', seq=benchmark_samples_all)
    output: 'benchmark/summary.tsv'
    shell:
        '''
        cat {input} > {output}
        '''

rule compile_benchmark:
    input: 'benchmark/{seq}/metaphlan2.tsv'
    output: 'benchmark/{seq,[.\w]+}.tsv'
    shell:
        r'''
        cd $(dirname {input})
        for i in *.tsv; do
          CLS=$(echo $i | rev | cut -d. -f2- | rev | cut -d. -f1)
          DB=$(echo $i | rev | cut -d. -f2- | rev | cut -d. -f2 -s)
          awk -v SEQ={wildcards.seq} -v DB="$DB" -v CLS="$CLS" 'BEGIN {{FS=OFS="\t"}} {{print SEQ, CLS, DB, $0}}' $i | tail -n +2;
        done > ../{wildcards.seq}.tsv
        '''

METHODS=['bracken',
'centrifuge',
'clark',
'clark_s',
'diamond',
'kaiju',
'kraken',
'kraken2',
'krakenhll',
'kslam',
'mmseqs2',
'megablast',
'prophyle',
'taxmaps',
]
rule compile_db_benchmark_all:
    input: expand('benchmark/db/{method}/refseqc.tsv', method=METHODS)
    output: 'benchmark/db/summary.tsv'
    shell:
        r'''
        for i in {input}; do
          awk -v CLS=$(basename $(dirname $i)) 'BEGIN {{FS=OFS="\t"}} {{print CLS, $0}}' $i | tail -n +2;
        done > {output}
        '''

rule compile_database_sizes:
    input: blast='db/refseqc/blast/genomic.nal'
    output: 'summary/database_sizes.tsv'
    shell:
        '''
        OUTPUT=$(readlink -f {output})
        pushd db/refseqc
        cd blast
        echo -e "megablast\trefseqc\t$(du -c genomic* taxdb* | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../bracken
        echo -e "bracken\trefseqc\t$(du -c KMER_DISTR.150* database150* | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../centrifuge
        echo -e "centrifuge\trefseqc\t$(du -c refseqc.*.cf | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../clark
        echo -e "clark\trefseqc\t$(du -c custom_0/db_central* | tail -n 1 | cut -f1)" >> "$OUTPUT"
        echo -e "clark_s\trefseqc\t$(du -c custom_0/T** | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../diamond
        echo -e "diamond\trefseqc\t$(du -c * | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../kaiju
        echo -e "kaiju\trefseqc\t$(du -c * | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../karp
        echo -e "karp\trefseqc\t$(du -c karp.tax karp.fasta rrna_karp.index | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../kraken
        echo -e "kraken\trefseqc\t$(du -c database.idx database.kdb taxonomy | tail -n 1 | cut -f1)" >> "$OUTPUT"

        cd ../kraken2
        echo -e "kraken2\trefseqc\t$(du -c *.k2d | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../krakenhll
        echo -e "krakenhll\trefseqc\t$(du -c database.idx database.kdb taxDB | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../kslam
        echo -e "kslam\trefseqc\t$(du -c database taxDB | tail -n 1 | cut -f1)" >> "$OUTPUT"
        # cd ../metaothello
        # echo -e "metaothello\trefseqc\t$(du -c id2tax names output.index | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../mmseqs2
        echo -e "mmseqs2\trefseqc\t$(du -c * | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../prophyle
        echo -e "prophyle\trefseqc\t$(du -c abv | tail -n 1 | cut -f1)" >> "$OUTPUT"
        cd ../taxmaps
        echo -e "taxmaps\trefseqc\t$(du -c refseqc.gem refseqc.len taxonomy.tbl | tail -n 1 | cut -f1)" >> "$OUTPUT"
        popd

        pushd db/gottcha
        echo -e "gottcha\tdefault\t$(du -c GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species.* GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.species.* *.dmp | tail -n 1 | cut -f1)" >> "$OUTPUT"
        popd

        pushd db/metaphlan2
        echo -e "metaphlan2\tdefault\t$(du -c db_v20 | tail -n 1 | cut -f1)" >> "$OUTPUT"
        popd

        pushd ~/efs/miniconda3/envs/metax/share/motus-2.0.1
        echo -e "motus2\tdefault\t$(du -c db_mOTU | tail -n 1 | cut -f1)" >> "$OUTPUT"
        popd

        pushd db/pathseq
        echo -e "pathseq\tdefault\t$(du -c pathseq_microbe* pathseq_taxonomy.db | tail -n 1 | cut -f1)" >> "$OUTPUT"
        popd

        pushd db/metaothello
        echo -e "metaothello\tdefault\t$(du -c 20161108 | tail -n 1 | cut -f1)" >> "$OUTPUT"
        popd
        '''

rule compile_classified_counts:
    input: expand('info/{sample}.metaphlan2.classified_count.txt', sample=samples_all),
           expand('info/{sample}.pathseq.classified_count.txt', sample=samples_all),
           expand('info/{sample}.kslam.{db}.classified_count.txt', sample=samples_all, db=['default', 'refseqc'])
    output: 'summary/classified_counts.tsv'
    shell:
        '''
        for i in {input}; do
           NAME="$(echo $(basename $i) | sed 's/.txt//' )"
           echo -e "$NAME\t$(cat $i)" >> {output}
        done
        '''

rule compile_reports:
    input: classified_counts='summary/classified_counts.tsv',
           total_reads='plotting/read_counts.tsv',
           classifer_outs=ALL_CLASSIFIERS_ALL
    output: reports='plotting/compiled_reports.tsv',
            ranksum='plotting/rank_abundances.tsv'
    shell:
        '''
        metax compile --tax-dir {config[TAXONOMY_DB]} --classified-counts {input.classified_counts} --read-counts {input.total_reads} --dir reports --config config/compile_config.yaml --rank-abundances {output.ranksum} {output.reports}
        '''
