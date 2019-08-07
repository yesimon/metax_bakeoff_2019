import argparse
import collections
import csv
import logging
import os
import re
import sys
import textwrap

import ncbitax
from metax import ioutil

log = logging.getLogger()

bad_prefixes = [
    'UnAmbiguouslyMapped_ds.frankengenome.mix',
    ]

prefixes = [
    'UnAmbiguouslyMapped_ds.7',
    'UnAmbiguouslyMapped_ds.buccal',
    'UnAmbiguouslyMapped_ds.cityparks',
    'UnAmbiguouslyMapped_ds.frankengenome',
    'UnAmbiguouslyMapped_ds.gut',
    'UnAmbiguouslyMapped_ds.hous1',
    'UnAmbiguouslyMapped_ds.hous2',
    'UnAmbiguouslyMapped_ds.nycsm',
    'UnAmbiguouslyMapped_ds.soil',
    'atcc_even.art',
    'atcc_staggered.art',
    'viral_genomes.art',
    'viral_refseqc',
    # 'viral_genomes.art',
]


ART_QNAME_ACCESSION = re.compile(r'^@?(.*)\.[1-9]\d*.*')


def match_fn_prefix(bname):
    for prefix in bad_prefixes:
        if bname.startswith(prefix):
            return False
    for prefix in prefixes:
        if bname.startswith(prefix):
            return prefix


class Processor:

    def __init__(self, output_dir, a2t):
        self.a2t = a2t
        self.output_dir = output_dir

    def try_process_file(self, fn):
        basename = os.path.basename(fn)
        prefix = match_fn_prefix(basename)
        if not prefix:
            return
        self.a2t.get_read_taxids(prefix)

        if re.match('.*\.kraken\.*', fn):
            self.process_kraken(fn)
        elif re.match('.*\.krakenhll\.*', fn):
            self.process_kraken(fn)
        elif re.match('.*\.diamond_tax\.*', fn):
            self.process_tsv(fn, sub='.tsv$')
        elif re.match('.*\.clark(_s)*\.*', fn):
            self.process_clark(fn)
        elif re.match('.*\.kslam\.*', fn):
            self.process_tsv(fn, sub='.tsv.gz$')
        elif re.match('.*\.metaothello\.*', fn):
            self.process_tsv(fn, sub='.tsv.bz2$')
        elif re.match('.*\.centrifuge\.*', fn):
            self.process_tsv(fn, taxid_col=2, sub='.txt.bz2$', skip_lines=1)
        elif re.match('.*\.kaiju\.*', fn):
            self.process_tsv(fn, qname_col=1, taxid_col=2, sub='.out$')
        elif re.match('.*\.taxmaps\.*', fn):
            self.process_tsv(fn, taxid_col=5, sub='.lca.')

    def output_path(self, bname):
        if not bname.endswith('.tsv'):
            bname = '{}.tsv'.format(bname)
        out_fn = os.path.join(self.output_dir, bname)
        return out_fn

    def output_newer(self, fn, bname):
        out_fn = self.output_path(bname)
        if not os.path.isfile(out_fn):
            return False
        if os.stat(out_fn).st_mtime > os.stat(fn).st_mtime:
            return True
        else:
            return False

    def print_hits(self, hits, bname):
        out_fn = self.output_path(bname)
        with open(out_fn, 'wt') as out_f:
            print('taxid', 'truth_taxid', 'count', sep='\t', file=out_f)
            for k, n in hits.most_common():
                taxid, truth_taxid = k
                print(taxid, truth_taxid, n, sep='\t', file=out_f)

    def process_kraken(self, fn):
        bname = os.path.basename(fn)
        bname = re.sub('.reads.gz$', '', bname)
        if self.output_newer(fn, bname):
            return

        hits = collections.Counter()
        with ioutil.compressed_open(fn, 'rt') as f:
            for line in f:
                parts = line.split('\t')
                classified = parts[0]
                qname = parts[1]
                taxid = int(parts[2])
                mo = ART_QNAME_ACCESSION.match(qname)
                if not mo:
                    logging.warning('No match for qname: %s', qname)
                accession = mo.group(1)
                truth_taxid = self.a2t.taxids[accession]
                hits[(taxid, truth_taxid)] += 1


        self.print_hits(hits, bname)

    def process_tsv(self, fn, qname_col=0, taxid_col=1, sub=None, skip_lines=0):
        bname = os.path.basename(fn)
        if sub:
            bname = re.sub(sub, '', bname)
        if self.output_newer(fn, bname):
            return

        hits = collections.Counter()
        with ioutil.compressed_open(fn, 'rt') as f:
            if skip_lines > 0:
                for _ in range(skip_lines):
                    next(f)
            for line in f:
                parts = line.split('\t')
                qname = parts[qname_col]
                # taxmaps
                try:
                    if parts[taxid_col] == '-':
                        taxid = 0
                    else:
                        taxid = int(parts[taxid_col])
                except:
                    print(line)
                    raise
                mo = ART_QNAME_ACCESSION.match(qname)
                if not mo:
                    logging.warning('No match for qname: %s', qname)
                try:
                    accession = mo.group(1)
                except:
                    print(line)
                    raise
                truth_taxid = self.a2t.taxids[accession]
                hits[(taxid, truth_taxid)] += 1

        self.print_hits(hits, bname)

    def process_clark(self, fn):
        bname = os.path.basename(fn)
        bname = re.sub('.csv$', '', bname)
        if self.output_newer(fn, bname):
            return

        hits = collections.Counter()
        with ioutil.compressed_open(fn, 'rt') as in_f:
            header = next(in_f)
            if '1st_assignment' in header:
                clark_s = True
                taxid_col = 3
            else:
                clark_s = False
                taxid_col = 2

            for row in csv.reader(in_f):
                qname = row[0]
                if row[taxid_col] == 'NA':
                    continue
                try:
                    taxid = int(row[taxid_col])
                except:
                    print(row)
                    raise
                mo = ART_QNAME_ACCESSION.match(qname)
                if not mo:
                    logging.warning('No match for qname: %s', qname)
                accession = mo.group(1)
                truth_taxid = self.a2t.taxids[accession]
                hits[(taxid, truth_taxid)] += 1

        self.print_hits(hits, bname)


class AccessionMap:
    def __init__(self, root_dir):
        self.root_dir = root_dir
        self.taxids = {}
        self.cached_fn_prefixes = set()

    def get_read_taxids(self, fn_prefix):
        if fn_prefix in self.cached_fn_prefixes:
            return
        map_fn = os.path.join(self.root_dir, '{}.tsv'.format(fn_prefix))
        with open(map_fn) as f:
            for line in f:
                parts = line.split('\t')
                if len(parts) < 2:
                    accession = line.rstrip()
                    self.taxids[accession] = None
                    log.warning('Taxid not found for accession: %s', accession)
                    continue
                base_accession, taxid = parts
                taxid = int(taxid)
                self.taxids[base_accession] = taxid
        self.cached_fn_prefixes.add(fn_prefix)


description = textwrap.dedent('''\
Calculate classification accuracy at the read level.
''')
def readclass_report(args):

    a2t = AccessionMap(os.path.join('fastq', 'accessions'))
    proc = Processor('readclass', a2t)
    files = []
    if args.files:
        files.extend(args.files)
    if args.dir:
        for path in os.listdir(args.dir):
            files.append(os.path.join(args.dir, path))
    for fn in files:
        proc.try_process_file(fn)



def add_command(subparsers):
    parser = subparsers.add_parser('readclass', description=description)
    parser.add_argument('--dir', help='Directory of files to process')
    parser.add_argument('--files', nargs='+', help='List of files to process')
    parser.add_argument('--taxids', help='Directory of taxids')
    # ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=readclass_report)
