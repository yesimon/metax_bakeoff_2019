"""
Create mmseqs taxonomy mapping from a fasta file

"""
import os
import sys
import time
import argparse
import logging

import ncbitax

from metax import ioutil

log = logging.getLogger()


def split_file(fp, marker):
    BLOCKSIZE = 4096
    current = ''
    for block in iter(lambda: fp.read(BLOCKSIZE), ''):
        current += block
        while 1:
            markerpos = current.find(marker)
            if markerpos == -1:
                break
            yield current[:markerpos]
            current = current[markerpos + len(marker):]
    yield current


def create_mmseqs_taxonomy(args):
    if args.action == 'accessions':
        for i, header in enumerate(split_file(args.infile, '\x00')):
            accessions = [x.split(' ')[0].split('.')[0] for x in header.split('\x01')]
            for a in accessions:
                print(a)
        sys.exit()

    tax_db = ncbitax.TaxonomyDb(tax_dir=args.tax_db, load_names=True, load_nodes=True)
    accession_lookup = {}
    start_time = time.time()
    log.info('Loading accession to taxid mapping')
    if args.genomic_joined:
        with ioutil.compressed_open(args.genomic_joined, 'rt') as f:
            for line in f:
                parts = line.split('\t')
                accession = parts[0]
                taxid = parts[1]
                taxid = int(taxid)
                accession_lookup[base_accession] = taxid
    elif args.accession2taxid:
        for fn in args.accession2taxid:
            with ioutil.compressed_open(fn, 'rt') as f:
                for line in f:
                    parts = line.split('\t')
                    base_accession = parts[0]
                    if base_accession == 'accession':
                        continue
                    try:
                        taxid = parts[2]
                    except:
                        print(line)
                    taxid = int(taxid)
                    accession_lookup[base_accession] = taxid
    log.info('Loaded accession to taxid mapping: %.2ss', time.time() - start_time)

    for i, header in enumerate(split_file(args.infile, '\x00')):
        accessions = [x.split(' ')[0].split('.')[0] for x in header.split('\x01')]
        taxids = []
        for a in accessions:
            try:
                taxid = accession_lookup[a]
            except KeyError:
                log.info('Taxid for accession not found: %s', a)
            taxids.append(taxid)
        lca = tax_db.coverage_lca(taxids)
        if not lca:
            log.error('LCA not found for taxids %s', taxids)
            lca = 1
        print("{}\t{}".format(i+1, lca), file=args.output)


def add_command(subparsers):
    parser = subparsers.add_parser('create-mmseqs-taxonomy', description='Create mmseqs taxonomy tsv from db_h.')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--output', nargs='?', type=argparse.FileType('wt'),
                        default=sys.stdout)
    parser.add_argument('--tax-db')
    parser.add_argument('--action', default='tsv')
    parser.add_argument('--genomic-joined', help='tsv of accession, taxid, gi')
    parser.add_argument('--accession2taxid', nargs='+', help='tsv of accession, taxid, gi')
    logging.basicConfig(level=logging.INFO)
    parser.set_defaults(func=create_mmseqs_taxonomy)
