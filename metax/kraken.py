import re

ranks = {
    'S': -1,
    'G': -2,
    'F': -3,
    'O': -4,
    'P': -5
    }


class KrakenTaxon:
    __slots__ = ('percent', 'cum_reads', 'unique_reads', 'rank', 'taxid', 'name', 'indent', 'parent')
    def __init__(self, percent, cum_reads, unique_reads, rank, taxid, name, indent):
        self.percent = percent
        self.cum_reads = cum_reads
        self.unique_reads = unique_reads
        self.rank = rank
        self.taxid = taxid
        self.name = name
        self.indent = indent
        self.parent = None

    def __repr__(self):
        return 'KrakenTaxon(percent={}, cum_reads={}, unique_reads={}, rank={}, taxid={}, name={})'.format(
            self.percent, self.cum_reads, self.unique_reads, self.rank, self.taxid, self.name)

    @classmethod
    def create_tree(cls, taxons):
        indent_stack = {}
        for t in taxons:
            if t.name == 'unclassified':
                yield t
                continue
            elif t.name == 'root':
                indent_stack[t.indent] = t
                yield t
            else:
                t.parent = indent_stack[t.indent - 1]
                indent_stack[t.indent] = t
                yield t


def read_kraken_report(f):
    for line in f:
        parts = line.strip().split('\t')
        percent = float(parts[0])
        cum_reads = int(parts[1])
        unique_reads = int(parts[2])
        rank = parts[3]
        taxid = int(parts[4])
        mo = re.search('^( *)(.*)', parts[5])
        indent = len(mo.group(1)) / 2
        yield KrakenTaxon(percent, cum_reads, unique_reads, rank, taxid, mo.group(2), indent)


def generate_tree(kraken_taxons):
    tree = KrakenTaxon.create_tree(kraken_taxons)
    return tree


def process_kraken_report(kraken_report_fn, taxid=None):
    with open(kraken_report_fn) as f:
        kraken_taxons = read_kraken_report(f)
        tree = list(generate_tree(kraken_taxons))

        root_taxon = tree[1]
        if taxid:
            genus_taxon = None
            species_taxon = None
            for t in tree:
                if taxid == t.taxid:
                    assert t.rank == 'S'
                    species_taxon = t

                    tp = t
                    while True:
                        tp = tp.parent
                        if tp.rank == 'G':
                            genus_taxon = tp
                        if tp.taxid == 1:
                            break
            percents = []
            if species_taxon is not None:
                species_percent = species_taxon.percent
            else:
                species_percent = 0
            if genus_taxon is not None:
                genus_percent = genus_taxon.percent
            else:
                genus_percent = 0
            try:
                percents = ['{:.4f}'.format(x/100.) for x in [species_percent, genus_percent, root_taxon.percent]]
            except:
                print(kraken_report_fn, taxid)
                raise
            return [str(taxid)] + percents
