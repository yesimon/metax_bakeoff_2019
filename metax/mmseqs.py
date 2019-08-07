import argparse
import random
import collections

import ncbitax
from metax import ioutil


def choose_best_random(counter):
    if len(counter) > 1:
        most_common = counter.most_common()
        n_best_hit = most_common[0][1]
        hits = [x for x in most_common if x[1] == n_best_hit]
        return random.choice(hits)
    elif len(counter) == 1:
        return counter.popitem()


def mmseqs_report(args):
    db = ncbitax.TaxonomyDb.from_args(args, load_nodes=True, load_names=True, load_merged=True)

    total_species_counter = collections.Counter()
    total_genus_counter = collections.Counter()

    taxid_map = {}

    last_read_id = None
    species_counter = collections.Counter()
    genus_counter = collections.Counter()
    with ioutil.compressed_open(args.input, 'rt') as in_f:
        for line in in_f:
            parts = line.rstrip().split('\t')
            read_id, taxid, _, rank, tax_name = parts
            taxid = int(taxid)
            if read_id != last_read_id:
                best = choose_best_random(species_counter)
                if best:
                    total_species_counter[best[0]] += 1
                # else:
                #     total_species_counter['unclassified'] += 1

                best = choose_best_random(genus_counter)
                if best:
                    total_genus_counter[best[0]] += 1
                # else:
                #     total_genus_counter['unclassified'] += 1
                species_counter = collections.Counter()
                genus_counter = collections.Counter()
            if rank == 'species':
                species_counter[taxid] += 1
            elif rank == 'genus':
                genus_counter[taxid] += 1

            taxid_map[taxid] = tax_name
            last_read_id = read_id

    total_reads = args.total_reads
    species_total = sum(x[1] for x in total_species_counter.items())
    genus_total = sum(x[1] for x in total_genus_counter.items())

    with ioutil.compressed_open(args.output, 'wt') as out_f:


        print('\t'.join(['0', 'unclassified', str((total_reads - species_total) / total_reads), 'species']), file=out_f)
        for taxid, n_reads in total_species_counter.most_common():
            abundance = n_reads / total_reads
            print('\t'.join(str(x) for x in [taxid, taxid_map[taxid], abundance, 'species']), file=out_f)

        print('\t'.join(['0', 'unclassified', str((total_reads - genus_total) / total_reads), 'genus']), file=out_f)
        for taxid, n_reads in total_genus_counter.most_common():
            abundance = n_reads / total_reads
            print('\t'.join(str(x) for x in [taxid, taxid_map[taxid], abundance, 'genus']), file=out_f)


def add_command(subparsers):
    parser = subparsers.add_parser('mmseqs-report')
    parser.add_argument('input', help='Input mmseqs taxonomy tsv output')
    parser.add_argument('output', help='Output report')
    parser.add_argument('--total-reads', required=True, type=int)
    ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=mmseqs_report)
