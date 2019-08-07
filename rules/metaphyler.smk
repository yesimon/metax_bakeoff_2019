METAPHYLER = config.get('METAPHYLER', 'metaphyler.pl')


rule metaphyler_single:
    input: 'fasta/{seq}.both.fasta'
    output: report='reports/{seq}.metaphyler.txt', reads='data/{seq}.metaphyler.tab'
    shadow: 'shallow'
    log: 'log/metaphyler/{seq}.log'
    shell:
        '''
        {METAPHYLER} {input} mp
        mv mp.classify.tab {output.reads}
        mv mp.class.tab {output.report}
        '''
