import argparse
import collections
import itertools
import textwrap
import logging

import ncbitax
from metax import ioutil


log = logging.getLogger()


def qname_taxids(db, in_f, out_f, taxid_ind):
    for qname, pairs in itertools.groupby(group_qname(in_f, taxid_ind), lambda x: x[0]):
        taxids = [x[1] for x in pairs]
        if 0 in taxids:
            lca = 0
        else:
            lca = db.coverage_lca(taxids)
            if not lca:
                lca = 0
        yield qname, lca


def group_qname(in_f, taxid_ind):
    for line in in_f:
        parts = line.split('\t')
        qname = parts[0]
        taxid = int(parts[taxid_ind])
        yield qname, taxid


def taxid_report(args):
    taxid_ind = args.taxid_column - 1

    db = ncbitax.TaxonomyDb.from_args(args, load_nodes=True, load_names=True, load_merged=True)
    in_f = ioutil.compressed_open(args.input)

    with in_f:
        with ioutil.compressed_open(args.output, 'wt') as out_f:
            lcas = collections.Counter(x[1] for x in qname_taxids(db, in_f, out_f, taxid_ind))
            for line in db.kraken_dfs_report(lcas):
                print(line, file=out_f)


description = textwrap.dedent('''\
Generate a kraken report for a tsv of read name, taxid
''')


def add_command(subparsers):
    parser = subparsers.add_parser('taxid-report', description=description)
    parser.add_argument('input')
    parser.add_argument('--output')
    parser.add_argument('--taxid-column', type=int, default=2)
    ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=taxid_report)
