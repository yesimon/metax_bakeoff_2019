import argparse

import ncbitax

def main():
    parser = argparse.ArgumentParser(description='Get all viral taxons')

    parser.add_argument('taxid_file')
    parser.add_argument('--tax-dir', required=False)
    args = parser.parse_args()

    db = ncbitax.TaxonomyDb(tax_dir=args.tax_dir, load_nodes=True, load_names=True, load_merged=True)

    with open(args.taxid_file) as f:
        for line in f:
            taxid = int(line)
            if db.is_viral(taxid):
                print(taxid)

if __name__ == '__main__':
    main()
