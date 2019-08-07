import argparse
import os
import sys
import collections

import ncbitax
from metax import ioutil


def metaothello_report(args):
    db = ncbitax.TaxonomyDb.from_args(args, load_nodes=True, load_names=True, load_merged=True)

    species_counter = collections.Counter()
    genus_counter = collections.Counter()
    with ioutil.compressed_open(args.input, 'rt') as in_f:
        for line in in_f:
            parts = line.split('\t')
            species_taxid, genus_taxid = int(parts[1]), int(parts[2])

            # Since metaothello db is old, some taxids get moved to subspecies, so let's attempt to fix
            if species_taxid > 1 and db.ranks.get(species_taxid) != 'species':
                path = db.parent_path(species_taxid)
                if path:
                    for p in path:
                        if db.ranks[p] == 'species':
                            species_taxid = p

            if genus_taxid > 1 and db.ranks.get(genus_taxid) != 'genus':
                path = db.parent_path(genus_taxid)
                if path:
                    for p in path:
                        if db.ranks[p] == 'genus':
                            genus_taxid = p

            species_counter[species_taxid] += 1
            genus_counter[genus_taxid] += 1

    with ioutil.compressed_open(args.output, 'wt') as out_f:
        neg_taxids = [x for x in set(species_counter.keys()) if x < 1]
        for taxid in neg_taxids:
            species_counter[0] += species_counter[taxid]
            del species_counter[taxid]
        species_total = sum(species_counter.values())
        for taxid, count in species_counter.items():
            if taxid == 0:
                name = 'unclassified'
            else:
                name = db.names[taxid]
            abundance = count / species_total
            print('\t'.join(str(x) for x in [taxid, name, abundance, 'species']), file=out_f)

        neg_taxids = [x for x in set(genus_counter.keys()) if x < 1]
        for taxid in neg_taxids:
            genus_counter[0] += genus_counter[taxid]
            del genus_counter[taxid]
        genus_total = sum(genus_counter.values())
        for taxid, count in genus_counter.items():
            if taxid == 0:
                name = 'unclassified'
            else:
                name = db.names[taxid]
            abundance = count / genus_total
            print('\t'.join(str(x) for x in [taxid, name, abundance, 'genus']), file=out_f)


def add_command(subparsers):
    parser = subparsers.add_parser('metaothello-report')
    parser.add_argument('input', help='Input metaothello read assignments')
    parser.add_argument('output', help='Output report')
    ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=metaothello_report)
