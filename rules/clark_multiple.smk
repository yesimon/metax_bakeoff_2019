import itertools

CLARK_PE_INPUT=list('fastq/%s.fastq' % ''.join(x) for x in itertools.product(samples_pe, paired_suffix))
CLARK_PE_OUTPUT=list('data/%s.fastq' % ''.join(x) for x in itertools.product(samples_pe, paired_suffix))

CLARK_SE_INPUT=list('fastq/%s.fastq' % x for x in samples_se)
CLARK_SE_OUTPUT=list('data/%s.fastq' % x for x in samples_se)

dbs = ['default_bv', 'refseqc']

samples_all = samples_se + samples_pe

rule clark_all:
    input: expand('reports/{sample}.{exe}.{db}.csv', exe=['clark', 'clark_s'], sample=samples_all, db=dbs)

rule clark_genus_all:
    input: expand('reports/{sample}.{exe}_genus.{db}.csv', exe=['clark', 'clark_s'], sample=samples_all, db=dbs)
    # input: expand('data/{sample}.{exe}_genus.{db}.csv', exe=['clark', 'clark_s'], sample=samples_all, db=dbs)

def clark_db_args(wildcards, output):
    if wildcards.db == 'default_bv':
        db = config.get('CLARK_BV_DB')
    elif wildcards.db == 'refseqc':
        db = config.get('CLARK_REFSEQC_DB')
    else:
        raise Exception

    if '_genus' in output[0]:
        db = db + '_1'
    else:
        db = db + '_0'

    targets = join(db, 'targets.txt')
    return '-D {} -T {}'.format(db, targets)

def clark_exe(wildcards):
    if wildcards.exe == 'clark':
        return config.get('CLARK')
    elif wildcards.exe == 'clark_s':
        return config.get('CLARK_S')
    raise Exception


# Run PE after SE in the same rule to reuse the db in cache as much as possible.
rule clark:
    input: se=CLARK_SE_INPUT,
           pe=CLARK_PE_INPUT
    output: aux_input_se=temp('aux/{exe}{rank}_input.{db}.txt'),
            aux_input_pe=[temp('aux/{exe}{rank}_input.{db}.1.txt'), temp('aux/{exe}{rank}_input.{db}.2.txt')],
            aux_output=temp('aux/{exe}{rank,(_genus)?}_results.{db}.txt'),
            files_se=list('data/%s.{exe}{rank,(_genus)?}.{db}.csv' % x for x in samples_se),
            files_pe=list('data/%s.{exe}{rank,(_genus)?}.{db}.csv' % x for x in samples_pe)
    log: log='log/clark/{exe}{rank}.{db}.log',
         time='time/clark/{exe}{rank}.{db}.log'
    params: exe=clark_exe,
            db=clark_db_args,
            output_files_se=lambda wildcards, output: [i[:-4] for i in output.files_se],
            output_files_pe=lambda wildcards, output: [i[:-4] for i in output.files_pe]
    wildcard_constraints:
        exe="(clark)|(clark_s)"
    threads: ALL_CORES
    shell:
        '''
        : > {output.aux_input_se}
        FILES=({input.se})
        for (( i=0; i<${{#FILES[@]}} ; i+=1 )) ; do
            echo "${{FILES[i]}}" >> {output.aux_input_se}
        done

        : > {output.aux_output}
        for i in {params.output_files_se}; do
          echo "$i" >> {output.aux_output}
        done

        /usr/bin/time -v -o {log.time} \
        {params.exe} {params.db} -n {threads} -O {output.aux_input_se} -R {output.aux_output} 2>&1 | tee {log.log}

        : > {output.aux_input_pe[0]}
        : > {output.aux_input_pe[1]}
        FILES=({input.pe})
        for (( i=0; i<${{#FILES[@]}} ; i+=2 )) ; do
            echo "${{FILES[i]}}" >> {output.aux_input_pe[0]}
            echo "${{FILES[i+1]}}" >> {output.aux_input_pe[1]}
        done

        : > {output.aux_output}
        for i in {params.output_files_pe}; do
            echo "$i" >> {output.aux_output}
        done


        /usr/bin/time -av -o {log.time} \
        {params.exe} {params.db} -n {threads} -P {output.aux_input_pe[0]} {output.aux_input_pe[1]} -R {output.aux_output} 2>&1 | tee -a {log.log}
        '''

def clark_abundance_exe(wildcards):
    dirname = clark_exe(wildcards)
    return os.path.join(os.path.dirname(dirname), 'getAbundance')

def clark_abundance_db_args(wildcards, output):
    if wildcards.db == 'default_bv':
        db = config.get('CLARK_BV_DB')
    elif wildcards.db == 'refseqc':
        db = config.get('CLARK_REFSEQC_DB')
    else:
        raise Exception

    if '_genus' in output[0]:
        db = db + '_1'
    else:
        db = db + '_0'

    return '-D {}'.format(os.path.dirname(db))

rule clark_abundance:
    input: 'data/{seq}.{exe}{rank,(_genus)?}.{db}.csv'
    output: 'reports/{seq}.{exe}{rank,(_genus)?}.{db}.csv'
    params: exe=clark_abundance_exe,
            db=clark_abundance_db_args
    wildcard_constraints:
        exe="(clark)|(clark_s)"
    shell:
        '''
        {params.exe} -F {input} {params.db} > {output}
        '''
