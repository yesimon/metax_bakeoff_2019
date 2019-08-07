#!/usr/bin/env python3
import argparse
import sys
import os
from pathlib import Path
from Bio import SeqIO
import gzip


def add_command(subparsers):
    parser = subparsers.add_parser('extract-gbff-rrna', description='Extract rRNA entries from GBFF files.')
    parser.add_argument('-i', '--input-dir', help='Output directory')
    parser.add_argument('-o', '--output-dir', help='Output directory')
    parser.set_defaults(func=extract_gbff_rrna)


def extract_gbff_rrna(args):
    output_dir = Path(args.output_dir)
    for path in Path(args.input_dir).iterdir():
        if not path.name.endswith('.gbff.gz'):
            continue

        output_path = output_dir / (path.name[:-8] + '.fa')

        rrna_found = 0
        with gzip.open(str(path), 'rt') as f:
            with output_path.open('wt') as of:
                for record in SeqIO.parse(f, 'genbank'):
                    for feature in record.features:
                        if feature.type != 'rRNA':
                            continue
                        prods = feature.qualifiers.get('product')
                        if not prods:
                            continue
                        for prod in prods:

                            sub_record = record[feature.location.start:feature.location.end]
                            sub_record.id = '{}:{}-{}'.format(record.id, feature.location.start, feature.location.end)
                            SeqIO.write(sub_record, of, 'fasta')
                            print(prod)
                            rrna_found += 1
