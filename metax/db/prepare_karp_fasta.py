#!/usr/bin/env python3
"""
Create Karp database input fasta from fasta in format >taxid|(taxid)| ...

Ignoring parens

The Karp format requires:
>(fasta reference id)  Bacteria;Proteobacteria;Gamma... down the tree levels for sequence
"""
from os.path import join
import os
import sys
import time
import argparse

import logging
from Bio import SeqIO

import ncbitax

log = logging.getLogger()


def parse_taxid(s):
    parts = s.split('|')
    taxid = int(parts[1])
    return taxid

def full_name(parents, names, taxid):
    path = []
    while taxid != 1:
        if taxid not in names:
            log.warning('Taxonomy id %s not found in names db.', taxid)
            return None
        path.append(names[taxid])
        taxid = parents[taxid]
    return ';'.join(reversed(path))


def add_command(subparsers):
    parser = subparsers.add_parser('prepare-karp-fasta', description='Create Karp database input FASTA')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--output', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('--tax-db')
    parser.add_argument('--genomic-joined', help='tsv of accession, taxid, gi')
    logging.basicConfig(level=logging.INFO)
    parser.set_defaults(func=prepare_karp_fasta)


def prepare_karp_fasta(args):
    tax_db = ncbitax.TaxonomyDb(tax_dir=args.tax_db, load_names=True, load_nodes=True)

    if args.genomic_joined:
        accession_lookup = {}
        with open(args.genomic_joined, 'rt') as f:
            for line in f:
                parts = line.split('\t')
                accession, taxid, gi = parts
                taxid = int(taxid)
                gi = int(gi)
                accession_lookup[accession] = taxid

        for seq in SeqIO.parse(args.infile, 'fasta'):
            seq_id = seq.id
            seq_id = seq_id.split(':')[0]
            taxid = accession_lookup.get(seq_id)
            if not taxid:
                log.info('Sequence ID %s not found in taxonomy db', seq_id)
                continue

            karp_name = full_name(tax_db.parents, tax_db.names, taxid)
            if not karp_name:
                continue
            seq.description = '\t'.join([seq.id, karp_name])
            SeqIO.write(seq, args.output, 'fasta')

        return


    for seq in SeqIO.parse(args.infile, 'fasta'):
        taxid = parse_taxid(seq.description)
        if taxid == 0:
            continue
        # seq.id = str(taxid)
        # seq.description = '\t'.join([str(taxid), full_name(tax_db.parents, tax_db.names, taxid)])
        # seq.description = '\t'.join([str(taxid), seq.description])
        karp_name = full_name(tax_db.parents, tax_db.names, taxid)
        if not karp_name:
            continue
        seq.description = '\t'.join([seq.id, karp_name])

        SeqIO.write(seq, args.output, 'fasta')
