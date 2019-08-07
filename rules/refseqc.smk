from os.path import join

NCBI_GENOMES = '/mnt/metax/db/ncbi_genomes'

rule download_refseq_genomes:
    output: join(NCBI_GENOMES, 'done')
    shell:
      '''
      ncbi-genome-download -F fasta,protein-fasta,genbank archaea,bacteria,viral
      touch {output}
      '''


rule refseqc_fasta:
    input: join(NCBI_GENOMES, 'done')
    output: fna='db/refseqc/fasta/genomic.fna.gz',
            faa='db/refseqc/fasta/protein.faa.gz',
            gbff='db/refseqc/fasta/genomic.gbff'
    shell:
      '''
      CWD=$(pwd)
      : > {output}

      cd {NCBI_GENOMES}/refseq

      DOMAINS=(archaea bacteria viral)
      for DOMAIN in ${{DOMAINS[@]}}; do
        pushd $DOMAIN
        awk -F "\t" '$12=="Complete Genome" && $11=="latest" {{print $1}}' assembly_summary.txt > ${{DOMAIN}}_selection.txt

        for i in $(awk -F "\t" '$12=="Complete Genome" {{print $1}}' assembly_summary.txt); do
          cat $i/*_genomic.fna.gz >> $CWD/{output.fna}
          cat $i/*_protein.faa.gz >> $CWD/{output.faa} || true
          {PIGZ} -dc $i/*_genomic.gbff.gz >> $CWD/{output.gbff}
        done
        popd
      done
      '''

rule refseqc_protein_uniq:
    input: 'db/refseqc/fasta/protein.faa.gz'
    output: 'db/refseqc/fasta/protein.uniq.faa.gz'
    threads: ALL_CORES
    resources: mem=ALL_MEMORY
    shell:
      r'''
      export LC_ALL=C
      {PIGZ} -dc {input} | fasta_formatter -t | sort -S {resources.mem}G -k1,1 -s -u | \
      awk 'BEGIN {{FS=OFS="\t"}} {{printf(">%s\n%s\n", $1, $2)}}' | {PIGZ} -9 > {output}
      '''


rule refseqc_fasta_decompressed:
    input: 'db/refseqc/genomic.fna.gz'
    output: 'db/refseqc/genomic.fna'
    shell:
        '''
        {PIGZ} -dk {input}
        '''

rule all_prot_accession2taxid:
    input: prot=join(TAXONOMY_DB, 'accession2taxid', 'prot.accession2taxid.gz'),
           dead_prot=join(TAXONOMY_DB, 'accession2taxid', 'dead_prot.accession2taxid.gz')
    output: join(TAXONOMY_DB, 'accession2taxid', 'prot.accession2taxid.vsorted')
    threads: ALL_CORES
    resources: mem=ALL_MEMORY
    shell:
      '''
      export LC_ALL=C
      cat <({PIGZ} -dc {input.prot} | tail -n +2) <({PIGZ} -dc {input.dead_prot} | tail -n +2) \
      | cut -f2,3 | sort -S {resources.mem}G -k1,1 -t$'\t' > {output}
      '''

# For taxmaps
rule all_nucl_accession2taxid:
    input: nucl=join(TAXONOMY_DB, 'accession2taxid', 'nucl_gb.accession2taxid.gz'),
           dead_nucl=join(TAXONOMY_DB, 'accession2taxid', 'dead_nucl.accession2taxid.gz')
    output: join(TAXONOMY_DB, 'accession2taxid', 'all_nucl.accession2taxid.sorted')
    threads: ALL_CORES
    resources: mem=ALL_MEMORY
    shell:
      '''
      export LC_ALL=C
      cat <({PIGZ} -dc {input.nucl} | tail -n +2) <({PIGZ} -dc {input.dead_nucl} | tail -n +2) \
      | cut -f2,3,4 | sort -S {resources.mem}G -k1,1 -t$'\t' > {output}
      '''

rule refseqc_genomic_headers_sorted:
    input: 'db/refseqc/fasta/genomic.fna.gz'
    output: 'db/refseqc/fasta/genomic.headers.sorted'
    resources: mem=ALL_MEMORY
    shell:
        '''
        export LC_ALL=C
        zcat {input} | fasta_formatter -t | cut -f1 | sed 's/ /\t/' | sort -k1,1 -s -S {resources.mem}G -t $'\t' > {output}
        '''

# For KrakenHLL
rule refseqc_genomic_headers_map:
    input: headers='db/refseqc/fasta/genomic.headers.sorted',
           a2t=join(TAXONOMY_DB, 'accession2taxid', 'all_nucl.accession2taxid.sorted')
    output: 'db/refseqc/fasta/genomic.map'
    shell:
        '''
        export LC_ALL=C
        join -t$'\t' <(sed 's/ /\t/' {input.headers}) {input.a2t} | cut -f1,3,2 > {output}
        '''

rule refseqc_genomic_gifna:
    input: fna='db/refseqc/fasta/genomic.fna.gz',
           a2t=join(TAXONOMY_DB, 'accession2taxid', 'all_nucl.accession2taxid.sorted')
    output: 'db/refseqc/fasta/genomic.gi.fna'
    resources: mem=ALL_MEMORY
    shell:
        '''
        export LC_ALL=C
        join -t$'\t' {input.a2t} <({PIGZ} -dc {input.fna} | fasta_formatter -t | sed 's/ /\t/' | sort -k1,1 -S {resources.mem}G -t $'\t') | \
        awk 'BEGIN {{FS=OFS="\t"}} {{printf(">gi|%s|%s %s\n%s\n", $3, $1, $4, $5)}} > {output}
        '''
