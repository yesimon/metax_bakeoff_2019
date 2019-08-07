#!/usr/bin/env python3
import argparse
import sys
import collections
import ncbitax


def create_metaothello_taxinfo(args):
    taxids = set()
    with open(args.taxids, 'rt') as f:
        for line in f:
            taxid = int(line)
            # accession, taxid, gi = line.split('\t')
            # taxid = int(taxid)
            # accession, taxid, gi = line.split('\t')
            # taxid = int(taxid)
            # gi = int(gi)
            taxids.add(taxid)

    tax_db = ncbitax.TaxonomyDb(tax_dir=args.tax_dir, load_nodes=True, load_names=True, load_merged=True)

    taxid_info = {}
    for taxid in taxids:
        inf = taxid_info[taxid] = {'species_taxid': taxid,
                                   'species_name': tax_db.names[taxid]}

        while True:
            taxid = tax_db.parents[taxid]
            if taxid == 1:
                break
            if tax_db.ranks[taxid] == 'genus':
                inf['genus_taxid'] = taxid
                inf['genus_name'] = tax_db.names[taxid]
                continue

            elif tax_db.ranks[taxid] == 'family':
                inf['family_taxid'] = taxid
                inf['family_name'] = tax_db.names[taxid]
                continue

            elif tax_db.ranks[taxid] == 'order':
                inf['order_taxid'] = taxid
                inf['order_name'] = tax_db.names[taxid]
                continue


            elif tax_db.ranks[taxid] == 'class':
                inf['class_taxid'] = taxid
                inf['class_name'] = tax_db.names[taxid]
                continue

            elif tax_db.ranks[taxid] == 'phylum':
                inf['phylum_taxid'] = taxid
                inf['phylum_name'] = tax_db.names[taxid]
                continue


            elif tax_db.ranks[taxid] == 'kingdom':
                inf['kingdom_taxid'] = taxid
                inf['kingdom_name'] = tax_db.names[taxid]
                continue

    def create_index(key):
        return {taxid: i for i, taxid in enumerate(inf.get(key + '_taxid', -1) for inf in taxid_info.values())}


    species_index = create_index('species')
    genus_index = create_index('genus')
    family_index = create_index('family')
    order_index = create_index('order')
    class_index = create_index('class')
    phylum_index = create_index('phylum')
    kingdom_index = create_index('kingdom')


    with open(args.output, 'wt') as f:
        header = '\t'.join(['Species_index','Species_ID','Species_name','Genus_index','Genus_ID','Genus_name','Family_index','Family_ID','Family_name','Order_index','Order_ID','Order_name','Class_index','Class_ID','Class_name','Phylum_index','Phylum_ID','Phylum_name'])

        print(header, file=f)

        for taxid in taxids:
            inf = taxid_info[taxid]
            sid = inf.get('species_taxid', -1)
            gid = inf.get('genus_taxid', -1)
            fid = inf.get('family_taxid', -1)
            oid = inf.get('order_taxid', -1)
            cid = inf.get('class_taxid', -1)
            pid = inf.get('phylum_taxid', -1)
            parts = [species_index[sid], sid, inf.get('species_name', 'NULL'),
                     genus_index[gid], gid, inf.get('genus_name', 'NULL'),
                     family_index[fid], fid, inf.get('family_name', 'NULL'),
                     order_index[oid], oid, inf.get('order_name', 'NULL'),
                     class_index[cid], cid, inf.get('class_name', 'NULL'),
                     phylum_index[pid], pid, inf.get('phylum_name', 'NULL'),
                     ]
            line = '\t'.join(str(x) for x in parts)
            print(line, file=f)


def add_command(subparsers):
    parser = subparsers.add_parser('create-metaothello-taxinfo', description='Create metaOthello taxinfo file.')
    parser.add_argument('taxids')
    parser.add_argument('--tax-dir')
    parser.add_argument('--output')
    parser.set_defaults(func=create_metaothello_taxinfo)
