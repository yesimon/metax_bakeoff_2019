import contextlib
import collections
import itertools
import tempfile
import operator

import ncbitax

from metax import ioutil

class BlastRecord(object):
    __slots__ = ('query_id', 'subject_id', 'percent_identity', 'aln_length', 'mismatch_count', 'gap_open_count', 'query_start',
                 'query_end', 'subject_start', 'subject_end', 'e_val', 'bit_score', 'gi', 'accession', 'taxids', 'scinames', 'comnames', 'title')

    def __init__(self, *args):
        self.query_id = None
        self.subject_id = None
        self.percent_identity = None
        self.aln_length = None
        self.mismatch_count = None
        self.gap_open_count = None
        self.query_start = None
        self.query_end = None
        self.subject_start = None
        self.subject_end = None
        self.e_val = None
        self.bit_score = None
        self.gi = None
        self.accession = None
        self.taxids = []
        self.scinames = []
        self.comnames = []
        self.title = []
        for attr, val in zip(self.__slots__, args):
            setattr(self, attr, val)


def blast_records(f):
    '''Yield blast m8 records line by line'''
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        for field in range(3, 10):
            parts[field] = int(parts[field])
        if len(parts) > 12:
            parts[12] = int(parts[12])
            parts[14] = [int(x) for x in parts[14].split(';')
                         if x != 'N/A']
        for field in (2, 10, 11):
            parts[field] = float(parts[field])
        # args = parts[:12]
        # extra = parts[12:]
        # args.append(extra)

        yield BlastRecord(*parts)

def process_blast_hits(db, hits, top_percent):
    '''Filter groups of blast hits and perform lca.

    Args:
      db: (TaxonomyDb) Taxonomy db.
      hits: []BlastRecord groups of hits.
      top_percent: (float) Only consider hits within this percent of top bit score.

    Return:
      (int) Tax id of LCA.
    '''
    hits = [hit for hit in hits if len(hit.taxids)]
    if len(hits) == 0:
        return

    best_score = max(hit.bit_score for hit in hits)
    cutoff_bit_score = (100 - top_percent) / 100 * best_score
    valid_hits = (hit for hit in hits if hit.bit_score >= cutoff_bit_score)
    valid_hits = list(valid_hits)
    # Sort requires realized list
    valid_hits.sort(key=operator.attrgetter('bit_score'), reverse=True)
    if valid_hits:
        tax_ids = list(itertools.chain(*(hit.taxids for hit in valid_hits)))
        return db.coverage_lca(tax_ids, lca_percent=None)


def blast_lca(db,
              m8_file,
              output=None,
              paired=False,
              min_bit_score=None,
              max_expected_value=None,
              top_percent=None):
    '''Calculate the LCA taxonomy id for groups of blast hits.

    Writes tsv output: query_id \t tax_id

    Args:
      db: (TaxonomyDb) Taxonomy db.
      m8_file: (io) Blast m8 file to read.
      output: (io) Output file.
      paired: (bool) Whether to count paired suffixes /1,/2 as one group.
      min_bit_score: (float) Minimum bit score or discard.
      max_expected_value: (float) Maximum e-val or discard.
      top_percent: (float) Only this percent within top hit are used.
    '''
    min_bit_score = min_bit_score if min_bit_score is not None else 50
    max_expected_value = max_expected_value if max_expected_value is not None else 0.01
    top_percent = top_percent if top_percent is not None else 10
    records = blast_records(m8_file)
    records = (r for r in records if r.e_val <= max_expected_value)
    records = (r for r in records if r.bit_score >= min_bit_score)
    if paired:
        records = (paired_query_id(rec) for rec in records)
    blast_groups = (v for k, v in itertools.groupby(records, operator.attrgetter('query_id')))
    hits = collections.Counter()
    for blast_group in blast_groups:
        blast_group = list(blast_group)
        tax_id = process_blast_hits(db, blast_group, top_percent)
        hits[tax_id] += 1
        if output:
            query_id = blast_group[0].query_id
            print(query_id, tax_id, sep='\t', file=output)
    return hits


def blast_report(args):

    tax_db = ncbitax.TaxonomyDb.from_args(args, load_nodes=True, load_names=True, load_merged=True)
    with contextlib.ExitStack() as ctx:
        if args.blast_report and not args.blast_lca:
            _, blast_lca_fn = tempfile.mkstemp('.blastn.lca.tsv')
        else:
            blast_lca_fn = args.blast_lca

        blast_m8_f = ctx.enter_context(ioutil.compressed_open(args.blast_m8, 'rt'))
        blast_lca_f = ctx.enter_context(ioutil.compressed_open(blast_lca_fn, 'wt'))
        hits = blast_lca(tax_db, blast_m8_f, blast_lca_f, min_bit_score=args.min_bit_score,
                         max_expected_value=args.max_expected_value, top_percent=args.top_percent)

        if not args.blast_report:
            return

        blast_report_f = ctx.enter_context(ioutil.compressed_open(args.blast_report, 'wt'))
        if blast_report_f:
            for line in tax_db.kraken_dfs_report(hits, total_reads=args.total_reads):
                print(line, file=blast_report_f)


def add_command(subparsers):
    parser = subparsers.add_parser('blast-report')

    parser.add_argument('blast_m8', help='Input BLAST m8 file.')
    parser.add_argument('--blast-lca', help='Per-read LCA result.')
    parser.add_argument('--blast-report', help='Kraken-like output report.')
    parser.add_argument('--total-reads', type=int, required=True, help='Total reads for generating unclassified fraction in the blast report')
    parser.add_argument('--min-bit-score', type=float)
    parser.add_argument('--max-expected-value', type=float)
    parser.add_argument('--top-percent', type=float)
    ncbitax.add_taxonomy_arguments(parser)
    parser.set_defaults(func=blast_report)
