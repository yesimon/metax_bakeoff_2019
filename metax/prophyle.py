import argparse
import collections
import itertools
import sys
import os
from os.path import join
import logging

import ncbitax

log = logging.getLogger()


def qname_taxids(db, in_sam_f, out_f):
    for qname, pairs in itertools.groupby(taxid_from_sam(in_sam_f), lambda x: x[0]):
        taxids = [x[1] for x in pairs]
        lca = db.coverage_lca(taxids)
        if not lca:
            lca = 0
        yield qname, lca


def taxid_from_sam(in_sam_f):
    for line in in_sam_f:
        if line.startswith('@'):
            # Skip headers
            continue
        parts = line.rstrip().split('\t')
        qname = parts[0]
        rname = parts[2]
        flag = int(parts[1])
        if flag & 0x4:
            continue

        try:
            taxid = int(rname.split('-')[1])
        except IndexError:
            print(line, file=sys.stderr)

        yield qname, taxid


def prophyle_report(args):
    total_reads = args.total_reads
    if args.paired:
        total_reads /= 2

    db = ncbitax.TaxonomyDb(tax_dir=args.tax_dir, load_nodes=True, load_names=True, load_merged=True)
    with open(args.input, 'rt') as in_sam_f:
        with open(args.output, 'wt') as out_f:
            lcas = collections.Counter(x[1] for x in qname_taxids(db, in_sam_f, out_f))

            for line in db.kraken_dfs_report(lcas, total_reads=int(total_reads)):
                print(line, file=out_f)


def add_command(subparsers):
    parser = subparsers.add_parser('prophyle-report')
    parser.add_argument('input')
    parser.add_argument('--output')
    parser.add_argument('--tax-dir', required=True)
    parser.add_argument('--total-reads', required=True, type=int)
    parser.add_argument('--paired', action='store_true')
    parser.set_defaults(func=prophyle_report)
