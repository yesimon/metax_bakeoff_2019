import itertools
import collections
import os

if kraken_execution == 'multiple':
    all_kraken_inputs = collections.defaultdict(list)
    all_kraken_reports = collections.defaultdict(list)
    all_kraken_reads = collections.defaultdict(list)
    all_kraken_pipes = collections.defaultdict(list)
    for db in dbs:
        kraken_samples = []
        for sample in samples:
            kraken_report = join('reports', "{sample}.kraken.{db}.txt".format(sample=sample, db=db))
            kraken_input = [join('fastq', "{sample}{pair}.fastq.bz2".format(sample=sample, pair=pair)) for pair in paired_suffix]
            try:
                itime = max(os.path.getmtime(kraken_input[0]), os.path.getmtime(kraken_input[1]))
                otime = os.path.getmtime(kraken_report)
            except OSError:
                kraken_samples.append(sample)
                continue
            if itime > otime:
                kraken_samples.append(sample)

        if not kraken_samples:
            kraken_samples = samples_pe.copy()

        all_kraken_inputs[db] = [
            join('fastq', "{sample}{pair}.fastq.bz2".format(pair=p, sample=s)) for s, p in itertools.product(kraken_samples, paired_suffix)]

        all_kraken_reports[db] = expand(
            join('reports', "{sample}.kraken.{db}.txt"), sample=kraken_samples, db=db)

        all_kraken_reads[db] = expand(
            join('data', "{sample}.kraken.{db}.reads.gz"), sample=kraken_samples, db=db)

        all_kraken_pipes[db] = expand(
            join('data', "{sample}.kraken.{db}.reads"), sample=kraken_samples, db=db)

kraken_multiple_s3_shell = """
trap rmpipe SIGINT SIGTERM EXIT

OUTPUT_PIPES=({params.pipes})

function rmpipe() {{
    trap - SIGTERM
    for (( i=0; i<${{#OUTPUT_PIPES[@]}} ; i+=1 )) ; do
      rm "${{OUTPUT_PIPES[i]}}"
    done
    pkill -P $$ || true
}}

for (( i=0; i<${{#OUTPUT_PIPES[@]}} ; i+=1 )) ; do
  mkfifo "${{OUTPUT_PIPES[i]}}"
done

{KRAKEN} --threads {threads} --taxonomy {params.taxonomy_nodes} --db-pipe --db-index <(aws s3 cp {params.idx} - | lz4 -d) --db-file <(aws s3 cp {params.kdb} - | lz4 -d) --index-size {params.index_size} --db-size {params.db_size} --paired --fastq-input --bzip2-compressed ${{OUTPUT_PIPES[@]/#/--output }} {input} | tee {log} &

OUTPUT_READS=({output.reads})
OUTPUT_REPORTS=({output.reports})
for (( i=0; i<${{#OUTPUT_PIPES[@]}} ; i+=1 )) ; do
cat "${{OUTPUT_PIPES[i]}}" | tee >(pigz --best > "${{OUTPUT_READS[i]}}") | {KRAKEN_FILTER} --taxonomy-nodes {params.taxonomy_nodes} --threshold {KRAKEN_FILTER_THRESHOLD} | {KRAKEN_REPORT} --taxonomy-nodes {params.taxonomy_nodes} --taxonomy-names {params.taxonomy_names} > "${{OUTPUT_REPORTS[i]}}"
done

wait
"""


rule kraken_full_multiple:
    input: all_kraken_inputs['full']
    output: reads=all_kraken_reads['full'],
            reports=all_kraken_reports['full']
    params: pipes=all_kraken_pipes['full'],
            idx=config['KRAKEN_FULL_IDX'],
            kdb=config['KRAKEN_FULL_KDB'],
            taxonomy_nodes=KRAKEN_TAXONOMY_NODES,
            taxonomy_names=KRAKEN_TAXONOMY_NAMES,
            index_size=8589934608,
            db_size=113890043500
    log: 'log/kraken/kraken.full.log'
    benchmark: 'benchmark/kraken/kraken.full.log'
    threads: ALL_CORES
    resources: mem=105
    shell:
        kraken_multiple_s3_shell

rule kraken_default_multiple:
    input: all_kraken_inputs['default']
    output: reads=all_kraken_reads['default'],
            reports=all_kraken_reports['default']
    params: pipes=all_kraken_pipes['default'],
            idx=config['KRAKEN_DEFAULT_IDX'],
            kdb=config['KRAKEN_DEFAULT_KDB'],
            taxonomy_nodes=KRAKEN_TAXONOMY_NODES,
            taxonomy_names=KRAKEN_TAXONOMY_NAMES,
            index_size=8589934608,
            db_size=102595705072
    log: 'log/kraken/kraken.default.log'
    benchmark: 'benchmark/kraken/kraken.default.log'
    threads: ALL_CORES
    resources: mem=190
    shell:
        kraken_multiple_s3_shell

rule kraken_mini_multiple:
    input: all_kraken_inputs['mini']
    output: reads=all_kraken_reads['mini'],
            reports=all_kraken_reports['mini']
    params: pipes=all_kraken_pipes['mini'],
            idx=config['KRAKEN_MINI_IDX'],
            kdb=config['KRAKEN_MINI_KDB'],
            taxonomy_nodes=KRAKEN_TAXONOMY_NODES,
            taxonomy_names=KRAKEN_TAXONOMY_NAMES,
            index_size=536870928,
            db_size=3597732208
    log: 'log/kraken/kraken.mini.log'
    benchmark: 'benchmark/kraken/kraken.mini.log'
    threads: ALL_CORES
    resources: mem=5
    shell:
        kraken_multiple_s3_shell
