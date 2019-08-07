VIROMESCAN = '/usr/local/var/metax/db/viromescan'

rule viromescan_all:
    input: expand('reports/{seq}.viromescan.tsv', seq=samples_pe)

rule viromescan:
    input: expand('fastq/{{seq}}{pair}.fastq.gz', pair=paired_suffix)
    output: report='reports/{seq}.viromescan.tsv'
    log: 'log/viromescan/{seq}.log'
    benchmark: 'benchmark/viromescan/{seq}.log'
    threads: 6
    shadow: 'shallow'
    resources: mem=10
    shell:
        '''
        env BMFILTER={VIROMESCAN}/tools/bmfilter SRPRISM={VIROMESCAN}/tools/srprism EXTRACT_FA={VIROMESCAN}/tools/extract_fullseq {VIROMESCAN}/viromescan.sh -1 {input[0]}  -2 {input[1]} -d virus_ALL -m /usr/local/var/metax/db -o viromescan -p {threads}
        mv viromescan/Species_level_results-Counts.txt {output.report}
        '''
