from os.path import join

BLAST_NT_SHM_DB = '/dev/shm/nt'

blast_rule = r'''\
export BLASTDB={params.db}
if [[ -s "{params.tmpout}" ]]; then

  # Just in case the query name contains a slash
  LAST_QNAME="$(python3 ~/efs/src/enc/blast-truncate.py {params.tmpout} | sed 's/\//\\\//g')"

  trap rmpipe SIGINT SIGTERM EXIT
  function rmpipe() {{
    trap - SIGTERM
    rm {params.pipe}
    pkill -P $$ || true
  }}

  mkfifo {params.pipe}

  cat {params.pipe} >> {params.tmpout} &

  fasta_formatter -t -i {input.fasta} | sed -n "/^$LAST_QNAME\t/,\$p" | awk '{{printf(">%s\n%s\n", $1, $2)}}' | \
    {params.blast} -db {params.db}/{params.db_prefix} -out {params.pipe} {params.max_hsps} -num_threads {threads} -task {params.task} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sgi sacc staxids sscinames scomnames stitle" 2> {log.err}
else
  {params.blast} -db {params.db}/{params.db_prefix} -query {input.fasta} -out {params.tmpout} {params.max_hsps} -num_threads {threads} -task {params.task} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sgi sacc staxids sscinames scomnames stitle" 2> {log.err}
fi

{PIGZ} -9 -p {threads} {params.tmpout}
mv {params.tmpout}.gz {output}
'''

blastn_nt_rsync_shm_rule = '''\
mkdir -p /dev/shm/nt
cd {params.db}
rsync -a {params.nt_set} /dev/shm/nt
rsync -a tax* /dev/shm/nt
rsync -a {params.nal} /dev/shm/nt/nt.nal
'''

blastdb_s3_shm_rule = '''\
cd /dev/shm
aws s3 cp {params.db} - | lbzip2 -dc | tar x
'''

BLAST_NT_S3_DB = config.get('BLAST_NT_S3_DB')
if BLAST_NT_S3_DB:
    def blast_db_s3(wildcards):
        if wildcards.db == 'refseqc':
            return config.get('BLAST_REFSEQC_S3_DB')
        elif wildcards.db == 'nt':
            return config.get('BLAST_NT_S3_DB')
        else:
            raise Exception

    rule blast_nt_s3_shm:
        output: touch('/dev/shm/{db}/{db}.done')
        params: db=blast_db_s3
        resources: mem=70
        priority: 78
        shell:
            blastdb_s3_shm_rule

else:
    BLAST_NT_DB = '/mnt/metax/db/blast/nt'
    rule blast_nt_rsync_shm:
        input: join(BLAST_NT_DB, 'nt.nal')
        output: touch('/dev/shm/nt/nt.done')
        params: db=BLAST_NT_DB,
                # nt_set='nt.0* nt.1* nt.2*',
                nt_set='nt.*',
                # nal='nt.00-30.nal'
                nal='nt.nal'
        resources: mem=70
        priority: 78
        shell:
            blastn_nt_rsync_shm_rule

rule makeblastdb_refseqc_nucl:
    input: fna='db/refseqc/fasta/genomic.fna',
           taxid_map='db/refseqc/fasta/genomic.joined'
    output: 'db/refseqc/blast/genomic.nal'
    params: out='db/refseqc/blast/genomic'
    benchmark: 'benchmark/db/megablast/refseqc.tsv'
    shadow: 'shallow'
    run:
        shell('dropcache')
        shell('makeblastdb -in {input.fna} -dbtype nucl -out {params.out} -taxid_map <(cut -f1-2 {input.taxid_map}) -parse_seqids',
              bench_record=bench_record)

def fasta_both_input(wildcards):
    if wildcards.seq in samples_pe:
        return 'fasta/{seq}.both.fasta'.format(seq=wildcards.seq)
    elif wildcards.seq in samples_se:
        return 'fasta/{}.fasta'.format(wildcards.seq)
    elif wildcards.seq in samples_fasta:
        return 'fasta/{}.fasta'.format(wildcards.seq)
    else:
        raise Exception

MEGABLAST_ALL = expand('reports/{sample}.megablast.{db}.txt', sample=samples_all, db=['nt'])
rule megablast_nt_all:
    input: MEGABLAST_ALL

# MEGABLAST_REFSEQC_ALL = expand('reports/{sample}.megablast.{db}.txt', sample=samples_all, db=['refseqc'])
MEGABLAST_REFSEQC_ALL = expand('reports/{sample}.megablast.{db}.txt', sample=samples_all, db=['refseqc'])
rule megablast_refseqc_all:
    input: MEGABLAST_REFSEQC_ALL


MEGABLAST_HG38_ALL = expand('data/hg38.{i}.megablast.{db}.tsv.gz', i=['{:02d}'.format(x) for x in range(16)], db=['refseqc'])
rule megablast_hg_all:
    input: MEGABLAST_HG38_ALL

MEGABLAST_HG38_NT_ALL = expand('data/hg38.{i}.megablast.{db}.tsv.gz', i=['{:02d}'.format(x) for x in range(32)], db=['nt'])
rule megablast_hg_nt_all:
    input: MEGABLAST_HG38_NT_ALL

BLAST_REFSEQC_DB = config.get('BLAST_REFSEQC_DB', '/dev/shm/refseqc')
BLAST_NT_DB = config.get('BLAST_NT_DB', '/dev/shm/nt')

def blast_db(wildcards):
    if wildcards.db == 'refseqc':
        return config.get('BLAST_REFSEQC_DB', '/dev/shm/refseqc')
    elif wildcards.db == 'nt':
        return config.get('BLAST_NT_DB', '/dev/shm/nt')
    else:
        raise Exception

rule megablast:
    input: fasta=fasta_both_input,
           db=ancient('/dev/shm/{db}/{db}.done')
    output: 'data/{seq}.megablast.{db,nt|refseqc}.tsv.gz'
    params: blast='blastn',
            db=blast_db,
            db_prefix='{db}',
            tmpout='tmp/{seq}.megablast.{db}.tmp',
            pipe='tmp/{seq}.megablast.{db}.pipe',
            task='megablast',
            max_hsps='-max_hsps 10'
    priority: 77
    resources: mem=8
    threads: 4
    log: log='log/megablast/{seq}.{db}.log',
         time='time/megablast/{seq}.{db}.log',
         err='log/megablast/{seq}.{db}.err'
    shell:
        blast_rule

# rule megablast_nt_benchmark:
#     input: fasta=fasta_both_input,
#             db=ancient('/dev/shm/nt/nt.done')
#     output: 'benchmark/data/{seq}.megablast.tsv.gz'
#     params: blast='blastn',
#             db=BLAST_NT_SHM_DB,
#             db_prefix='nt',
#             tmpout='tmp/{seq}.megablast.tmp',
#             pipe='tmp/{seq}.megablast.pipe',
#             task='megablast',
#             max_hsps='-max_hsps 10'
#     priority: 77
#     resources: mem=4
#     threads: ALL_CORES
#     log: log='benchmark/log/megablast/{seq}.log',
#          time='benchmark/time/megablast/{seq}.log'
#     benchmark: repeat('benchmark/{seq}/megablast.{db}.log', 2)
#     run:
#         if benchmark_i == 0:
#             shell('dropcache')
#         shell(blast_rule)

rule megablast_report:
    input: data='data/{seq}.megablast.{db}.tsv.gz',
           total_reads='info/{seq}.total_reads.txt'
    output: lca='data/{seq}.megablast.lca.{db}.tsv.gz',
            report='reports/{seq}.megablast.{db}.txt'
    params: db=config.get('TAXONOMY_DB')
    shell:
        '''
        metax blast-report --total-reads $(cat {input.total_reads}) --blast-lca {output.lca} --blast-report {output.report} --tax-dir {params.db} {input.data}
        '''
