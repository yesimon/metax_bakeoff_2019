fastq_pairing = config.get('fastq_pairing', 'interleaved' )

# Raw fastqs may be improperly paired - process and save only properly paired
# (id in both fastqs) reads.
rule proper_pair_fastq:
    input: expand('raw_fastq/{{seq}}{pair}.fastq.gz', pair=paired_suffix)
    output: expand('fastq/{{seq}}{pair}.fastq.gz', pair=paired_suffix)
    shell:
        r'''
        join -t$'\t' <(zcat {input[0]} | paste - - - - | sort -t$'\t' -k1,1) <(zcat {input[1]} | paste - - - - | sort -t$'\t' -k1,1) | \
        sort -t$'\t' -k1,1V  | \
        awk 'BEGIN {{FS=OFS="\t"}} {{ \
            printf("%s/1\n%s\n+\n%s\n", $1, $2, $4) | "pigz -9 > {output[0]}"; \
            printf("%s/2\n%s\n+\n%s\n", $1, $5, $7) | "pigz -9 > {output[1]}" \
        }}'
        '''


# Remove X nucleotides from Mavromatis*
rule clean_fasta_bz2:
    input: 'raw_fastq/{seq}.fasta.bz2'
    output: 'fastq/{seq}.fasta.bz2'
    shell:
        '''
        lbzip2 -dc {input} | sed -e '/^[^>]/s/X/N/g' | lbzip2 > {output}
        '''


rule pair_fastq_bz2:
  input: expand('fastq/{{sample}}{pair}.fastq.bz2', pair=paired_suffix)
  output: 'fastq/{sample}.both.fastq.bz2'
  # output: temp('fastq/{sample}.both.fastq.bz2')
  run:
    if fastq_pairing == 'interleaved':
        shell("""
        echo {input[0]}
        echo {input[1]}
        paste <(paste - - - - < <(lbzcat {input[0]})) \
              <(paste - - - - < <(lbzcat {input[1]})) \
              | tr '\t' '\n' \
              | lbzip2 > {output}
        """)
    elif fastq_pairing == 'interleaved_add_suffix':
        shell(r"""
        echo {input[0]}
        echo {input[1]}
        paste <(paste - - - - < <(awk '{{if (NR % 4 == 1) {{printf "%s/1\n", $1}} \
                                         else if (NR % 4 == 3){{print "+"}} else {{print $0}} }}' <(lbzcat {input[0]}))) \
              <(paste - - - - < <(awk '{{if (NR % 4 == 1) {{printf "%s/2\n", $1}} \
                                         else if (NR % 4 == 3){{print "+"}} else {{print $0}} }}' <(lbzcat {input[1]}))) \
              | tr '\t' '\n' \
              | lbzip2 > {output}
        """)
    elif fastq_pairing == 'cat':
        shell("""
        echo {input[0]}
        echo {input[1]}
        lbzcat {input[0]} {input[1]} | lbzip2 > {output}
        """)


rule unbzip2_fastq:
    input: 'fastq/{fastq}.fastq.bz2'
    output: 'fastq/{fastq}.fastq'
    # output: temp('fastq/{fastq}.fastq')
    shell:
        """
        echo {input}
        lbzip2 --decompress --force --keep {input[0]}
        # If original file is empty
        touch {output}
        """


rule unbzip2_fasta:
    input: 'fastq/{fastq}.fasta.bz2'
    output: 'fastq/{fastq}.fasta'
    # output: temp('fastq/{fastq}.fasta')
    shell:
      """
      echo {input}
      lbzip2 --decompress --force --keep {input[0]}
      # If original file is empty
      touch {output}
      """

rule fastqc:
    input: 'fastq/{sample}.both.fastq.gz'
    output: 'qc/{sample}_fastqc.html', 'qc/{sample}_fastqc.zip'
    shell:
        '''
        fastqc -o qc -f fastq <(zcat {input})
        '''

rule fastq_to_fasta:
    input: 'fastq/{seq}.fastq'
    # output: temp('fasta/{seq}.fasta')
    output: 'fasta/{seq}.fasta'
    shell:
        '''
        seqtk seq -a {input} > {output}
        '''

rule fasta_to_fastq:
    input: 'fastq/{seq}.fasta'
    output: 'fastq/{seq}.fastq'
    # output: temp('fastq/{seq}.fastq')
    shell:
      '''
      seqtk seq -F F {input} > {output}
      '''


rule cd_hit_dup_all:
    input: expand('fastq/cdhit/{sample}.1.fastq.clstr', sample=samples_fastq)


def cd_hit_dup_args(wildcards, input):
    if wildcards.seq in samples_pe:
        o = 'fastq/cdhit/{}.1.fastq'.format(wildcards.seq)
        o2 = 'fastq/cdhit/{}.2.fastq'.format(wildcards.seq)
        return ' -i {} -i2 {} -o {} -o2 {}'.format(input[0], input[1], o, o2)
    elif wildcards.seq in samples_se:
        o = 'fastq/cdhit/{}.1.fastq'.format(wildcards.seq)
        return ' -i {} -o {}'.format(input, o)
    else:
        raise Exception

rule cd_hit_dup:
    input: fastq_input
    output: 'fastq/cdhit/{seq}.1.fastq.clstr'
    params: ofastq1='fastq/cdhit/{seq}.1.fastq',
            ofastq2='fastq/cdhit/{seq}.2.fastq',
            args=cd_hit_dup_args,
            exe=config.get('CD_HIT_DUP', 'cd-hit-dup')
    log: log='log/cdhit/{seq}.log',
         time='time/cdhit/{seq}.log'
    threads: 1
    shell:
        '''
        /usr/bin/time -v -o {log.time} \
        {params.exe} {params.args} 2>&1 | tee {log}
        rm -f {params.ofastq1} {params.ofastq2}
        '''
