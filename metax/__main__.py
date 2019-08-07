#!/usr/bin/env python3
import argparse
import os
import csv
import collections
import re
import sys

import metax.blast_report
import metax.reporter
import metax.prophyle
import metax.metaothello
import metax.merge_jellyfish
import metax.mmseqs
import metax.taxid_report
import metax.readclass
import metax.db.prepare_karp_fasta
import metax.db.create_metaothello_taxinfo
import metax.db.extract_gbff_rrna
import metax.db.create_mmseqs_taxonomy
import metax.db.create_centrifuge_map
import metax.merge_sam


def main():
    parser = argparse.ArgumentParser(description='metax scripts and utilities')

    subparsers = parser.add_subparsers()
    metax.merge_sam.add_command(subparsers)
    metax.blast_report.add_command(subparsers)
    metax.reporter.add_command(subparsers)
    metax.metaothello.add_command(subparsers)
    metax.merge_jellyfish.add_command(subparsers)
    metax.prophyle.add_command(subparsers)
    metax.mmseqs.add_command(subparsers)
    metax.taxid_report.add_command(subparsers)
    metax.db.prepare_karp_fasta.add_command(subparsers)
    metax.db.create_metaothello_taxinfo.add_command(subparsers)
    metax.db.extract_gbff_rrna.add_command(subparsers)
    metax.db.create_mmseqs_taxonomy.add_command(subparsers)
    metax.db.create_centrifuge_map.add_command(subparsers)
    metax.readclass.add_command(subparsers)

    args = parser.parse_args()
    func = getattr(args, 'func', None)
    if func:
        args.func(args)
    else:
        print('No command selected')


if __name__ == '__main__':
    main()
