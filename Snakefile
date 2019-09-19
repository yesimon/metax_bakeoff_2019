shell.executable('/bin/bash')

import collections
import csv
import itertools
import os
from pathlib import Path
import re
import subprocess
import textwrap
from pandas.core.common import flatten
import psutil
ALL_MEMORY = int((psutil.virtual_memory().total / 1024 ** 3) - 1)
MEM_TMPDIR = config.get('MEM_TMPDIR', '/dev/shm')

configfile: 'config/config.yaml'
include: 'rules/config.smk'

TMPDIR = config.get('TMPDIR', '/tmp')
PIGZ = config.get('PIGZ', 'pigz')
SORT = config.get('SORT', 'sort')
PICARD = config.get('PICARD', 'picard')
SAMTOOLS = config.get('SAMTOOLS', 'samtools')
DROPCACHE = config.get('DROPCACHE', 'dropcache')

samples_se = [str(x) for x in set(config.get('samples_se', [])) - set(config.get('exclude_samples_pe', []))]
samples_pe = [str(x) for x in set(config.get('samples_pe', [])) - set(config.get('exclude_samples', []))]
samples_fastq = samples_se + samples_pe
samples_fasta = [str(x) for x in set(config.get('samples_fasta', []))]
samples_all = samples_fastq + samples_fasta

benchmark_samples_se = list(sorted(config.get('benchmark_samples_se')))
benchmark_samples_pe = list(sorted(config.get('benchmark_samples_pe')))
benchmark_samples_all = benchmark_samples_se + benchmark_samples_pe

paired_suffix = config.get('paired_suffix', ('_1', '_2'))

TAXONOMY_DB = config.get('TAXONOMY_DB')


def fastq_bz2_input(wildcards):
    if wildcards.seq in samples_pe:
        return expand('fastq/{seq}{pair}.fastq.bz2', seq=wildcards.seq, pair=paired_suffix)
    elif wildcards.seq in samples_se:
        return 'fastq/{}.fastq.bz2'.format(wildcards.seq)
    elif wildcards.seq in samples_fasta:
        return 'fastq/{}.fastq.bz2'.format(wildcards.seq)
    else:
        raise Exception

def fastx_bz2_input(wildcards):
    if wildcards.seq in samples_pe:
        return expand('fastq/{seq}{pair}.fastq.bz2', seq=wildcards.seq, pair=paired_suffix)
    elif wildcards.seq in samples_se:
        return 'fastq/{}.fastq.bz2'.format(wildcards.seq)
    elif wildcards.seq in samples_fasta:
        return 'fastq/{}.fasta.bz2'.format(wildcards.seq)
    else:
        raise Exception

def fastq_both_input(wildcards):
    if wildcards.seq in samples_pe:
        return 'fastq/{}.both.fastq'.format(wildcards.seq)
    elif wildcards.seq in samples_se:
        return 'fastq/{}.fastq'.format(wildcards.seq)
    elif wildcards.seq in samples_fasta:
        return 'fastq/{}.fastq'.format(wildcards.seq)
    else:
        raise Exception

def fastx_both_input(wildcards):
    if wildcards.seq in samples_pe:
        return 'fastq/{}.both.fastq'.format(wildcards.seq)
    elif wildcards.seq in samples_se:
        return 'fastq/{}.fastq'.format(wildcards.seq)
    elif wildcards.seq in samples_fasta:
        return 'fastq/{}.fasta'.format(wildcards.seq)
    else:
        raise Exception


def fastq_input(wildcards):
    if wildcards.seq in samples_pe:
        return expand('fastq/{seq}{{pair}}.fastq'.format(seq=wildcards.seq), pair=paired_suffix)
    elif wildcards.seq in samples_se:
        return 'fastq/{}.fastq'.format(wildcards.seq)
    elif wildcards.seq in samples_fasta:
        return 'fastq/{}.fastq'.format(wildcards.seq)
    else:
        raise Exception

def fastx_input(wildcards):
    if wildcards.seq in samples_pe:
        return expand('fastq/{seq}{{pair}}.fastq'.format(seq=wildcards.seq), pair=paired_suffix)
    elif wildcards.seq in samples_se:
        return 'fastq/{}.fastq'.format(wildcards.seq)
    elif wildcards.seq in samples_fasta:
        return 'fastq/{}.fasta'.format(wildcards.seq)
    else:
        raise Exception

include: 'rules/download.smk'
include: 'rules/bam_to_fastq.smk'
include: 'rules/prepare_fastq.smk'
include: 'rules/refseqc.smk'

include: 'rules/blast.smk'
include: 'rules/bracken.smk'
include: 'rules/centrifuge.smk'
clark_execution = config.get('CLARK_EXECUTION', 'single')
if clark_execution == 'multiple':
    include: 'rules/clark_multiple.smk'
elif clark_execution == 'single':
    include: 'rules/clark.smk'
include: 'rules/diamond.smk'
include: 'rules/acdiamond.smk'
include: 'rules/gottcha.smk'
include: 'rules/kaiju.smk'
include: 'rules/karp.smk'
include: 'rules/kslam.smk'
kraken_execution = config.get('KRAKEN_EXECUTION', 'single')
if kraken_execution == 'multiple':
    include: 'rules/kraken_multiple.smk'
elif kraken_execution == 'single':
    include: 'rules/kraken.smk'
include: 'rules/kraken2.smk'
include: 'rules/krakenhll.smk'
include: 'rules/metaothello.smk'
include: 'rules/metaphlan2.smk'
include: 'rules/mmseqs2.smk'
include: 'rules/motus.smk'
include: 'rules/pathseq.smk'
include: 'rules/prophyle.smk'
include: 'rules/taxmaps.smk'

ALL_CLASSIFIERS_ALL = list(flatten([
    MEGABLAST_ALL, CENTRIFUGE_ALL, CLARK_ALL, DIAMOND_ALL, GOTTCHA_ALL, KAIJU_ALL, KSLAM_ALL, KRAKEN_ALL, KRAKEN2_ALL, KRAKENHLL_ALL, METAOTHELLO_ALL, MMSEQS2_ALL, MOTUS_ALL, PATHSEQ_ALL, PROPHYLE_ALL, TAXMAPS_ALL]))
print(ALL_CLASSIFIERS_ALL)

include: 'rules/reports.smk'


rule all:
  input: MEGABLAST_ALL, CENTRIFUGE_ALL, CLARK_ALL, DIAMOND_ALL, GOTTCHA_ALL, KAIJU_ALL, KSLAM_ALL, KRAKEN_ALL, KRAKEN2_ALL, KRAKENHLL_ALL, METAOTHELLO_ALL, MMSEQS2_ALL, MOTUS_ALL, PATHSEQ_ALL, PROPHYLE_ALL, TAXMAPS_ALL
