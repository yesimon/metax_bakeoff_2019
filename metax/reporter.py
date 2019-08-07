#!/usr/bin/env python3
import argparse
import os
import csv
import logging
from pathlib import Path
import collections
import re
import sys

import yaml

import ncbitax

log = logging.getLogger()

def load_config(conf_path):
    if conf_path:
        with open(conf_path) as f:
            conf = yaml.safe_load(f)
    else:
        conf = {}
    return conf


KRAKEN_CODE_BACK = {
    'K': 'kingdom',
    'P': 'phylum',
    'C': 'class',
    'O': 'order',
    'F': 'family',
    'G': 'genus',
    'S': 'species',
    'D': 'superkingdom',
}


class ReportTableWriter:
    def __init__(self, outf):
        self.outf = outf
        fieldnames = ['sample', 'classifier', 'database', 'taxid', 'name', 'cum_abundance', 'cls_cum_abundance',
                      'unique_abundance', 'cls_unique_abundance', 'classrank', 'rank']
        self.writer = csv.DictWriter(outf, fieldnames=fieldnames, delimiter='\t',
                                     extrasaction='ignore')
        self.writer.writeheader()

    def write(self, d):
        self.writer.writerow(d)


class ReportProcessor:
    def __init__(self, db, conf, matchers=None):
        self.db = db
        self.matchers = matchers or Matchers.create(self)
        self.name_to_taxid = {}
        self.name_lower_to_taxid = {}
        self.name_lower_to_original = {}
        # self.scientific_name_to_taxid = {}
        self.species_to_taxid = {}
        self.conf = conf
        self.rank_abundances = collections.defaultdict(float)
        self.parent_path_cache = {}

        for taxid, names in self.db.names.items():
            # Don't let merged taxids overwrite the name -> taxid mapping
            if taxid in self.db.merged:
                continue
            for name in names:
                self.name_to_taxid[name] = taxid
                self.name_lower_to_taxid[name.lower()] = taxid
                self.name_lower_to_original[name.lower()] = name
        # for taxid, sciname in self.db.scientific_names.items():
        #     self.scientific_name_to_taxid[sciname] = taxid


        # for name, taxids in self.name_to_taxid.items():
        #     for taxid in taxids:
        #         if self.db.ranks[taxid] == 'species':
        #             self.name_to_taxid[name][0] = taxid
        self.name_failed = set(['classified_above_species'])


    def get_taxid(self, name, no_taxid=False):
        '''Set no_taxid to true if no taxid is expected for the name. Therefore don't output a warning.'''
        taxid = self.name_to_taxid.get(name)
        if taxid is not None:
            return taxid
        else:
            if name in self.name_failed:
                return
            self.name_failed.add(name)
            if not no_taxid:
                log.warning("%s: Couldn't find taxid for name '%s'", self.filename, name)

    def set_taxid(self, name, d):
        no_taxid = d.pop('no_taxid', False)
        res = self.get_taxid(name, no_taxid=no_taxid)
        if res:
            d['taxid'] = res

    def fn_prefix_info(self, fn_prefix):
        fn_prefix = fn_prefix.rstrip('.')
        info = self.conf['samples'].get(fn_prefix, {})
        if not info:
            log.debug('sample not found in config %s', fn_prefix)
        return info

        # if info:
        #     return info
        # else:
        #     log.error('sample not found in config %s', fn_prefix)

    def process_file(self, filename, in_f):
        filename = Path(filename)
        self.filename = filename
        tup = self.matchers.match(filename)
        success = tup[0]
        if not success:
            reason = tup[1]
            if reason == 'no_matcher':
                print('file failed - no matcher -  {}'.format(filename))
            return

        parser, mo = tup[1:]

        fn_prefix = mo.group('fn_prefix')
        info = self.fn_prefix_info(fn_prefix)

        # if info is None:
        #     return

        dbase = {'sample': info.get('sample', fn_prefix),
                 'paired': info.get('paired', False)}

        kwargs = {}
        try:
            database = mo.group('database')
            kwargs['database'] = database
        except IndexError:
            pass

        taxa = []
        try:
            for d in parser.parse(in_f, dbase, **kwargs):
                # Something failed in parsing, don't output taxa rows
                if not d:
                    return
                if not d['sample']:
                    print(in_f)

                classrank = getattr(parser, 'CLASSRANK', 'all')
                d.setdefault('classrank', classrank)
                if 'taxid' not in d and d['name'] not in parser.IGNORED_NAMES:
                    self.set_taxid(d['name'], d)

                if 'name' not in d and 'taxid' in d:
                    try:
                        d['name'] = self.db.names[d['taxid']][0]
                    except:
                        pass
                taxid = d.get('taxid')
                if taxid and not d.get('rank'):
                    d['rank'] = self.db.ranks.get(taxid)
                self.label_subspecies_ranks(d)
                taxa.append(d)

        except Exception as e:
            log.exception('Parse failed for file %s', in_f.name)

        return taxa

    def sum_rank_abundances(self, taxa):
        for d in taxa:
            try:
                if d.get('rank') in ['genus', 'species', 'subspecies']:
                    self.rank_abundances[(d['sample'], d['classifier'], d.get('database', 'default'), d['rank'])] += d['cum_abundance']
            except:
                print(d)
                raise

    def label_subspecies_ranks(self, d):
        # Label all children of species rank for total subspecies cum abundance and labelling 'below_species'
        taxid = d.get('taxid')
        if not taxid or taxid == 1:
            return
        path = self.db.parent_path(int(d['taxid']), cache=self.parent_path_cache, warn_missing=False)
        if not path:
            return
        for i, path_taxid in enumerate(path):
            if self.db.ranks[path_taxid] == 'species':
                if i == 0:
                    d['rank'] == 'subspecies'
                else:
                    # Sometimes past subspecies get moved in taxonomy tree to non-subspecies
                    if d['rank'] == 'subspecies':
                        log.info('Preexisting subspecies not direct child of species: %s', d['taxid'])
                    else:
                        d['rank'] == 'below_species'
                return

    def create_unclassified_taxon(self, d, unclassified_abundance=None):
        sample, classifier, db = d['sample'], d['classifier'], d['database']


        if unclassified_abundance is None:
            total_reads = self.total_reads.get(sample)
            if total_reads is None:
                log.warning('Missing reads count for sample: %s', sample)
                return

            classified_reads = self.classified_counts.get((sample, classifier, db))
            if classified_reads is None:
                log.warning('Missing classified count for file: %s', (sample, classifier, db))
                return

            unclassified_count = self.total_reads[sample] - self.classified_counts[(sample, classifier, db)]

            abundance = unclassified_count / self.total_reads[sample]
            if unclassified_count < 0:
                print(sample, classifier, db, self.total_reads[sample], unclassified_count)
        else:
            abundance = unclassified_abundance
            unclassified_count = None
        return {
            'sample': sample,
            'database': db,
            'classifier': classifier,
            'name': 'unclassified',
            'taxid': 0,
            'cum_abundance': abundance,
            'unique_count': unclassified_count,
            'classrank': 'all'
        }


class Matchers:
    def __init__(self):
        self.matchers = {}

    def add(self, matcher_name, matcher):
        self.matchers[matcher_name] = matcher

    def match(self, filepath):
        for matcher_name, parser in self.matchers.items():
            active = getattr(parser, 'active', True)
            success, mo = parser.match_filepath(filepath)
            if success:
                if not active:
                    return False, 'inactive'
                return True, parser, mo
            else:
                if mo == 'do_not_parse':
                    return False, 'do_not_parse'
        return False, 'no_matcher'

    @classmethod
    def create(cls, processor):
        m = cls()
        m.add('bracken', BrackenParser(processor))
        m.add('bracken_genus', BrackenGenusParser(processor))
        m.add('centrifuge', CentkrakenParser(processor, fix_species_group=True))
        m.add('centrifuge_raw', CentrifugeParser(processor))
        m.add('clark', ClarkParser(processor))
        m.add('clark_s', ClarkSParser(processor))
        m.add('clark_genus', ClarkGenusParser(processor))
        m.add('clark_s_genus', ClarkSGenusParser(processor))
        m.add('diamond', DiamondKrakenParser(processor))
        m.add('kaiju', KaijuParser(processor))
        m.add('kaiju_genus', KaijuGenusParser(processor))
        m.add('kraken', KrakenParser(processor))
        m.add('kraken2', Kraken2Parser(processor))
        m.add('krakenhll', KrakenhllParser(processor))
        m.add('gottcha', GottchaParser(processor))
        m.add('gottcha_genus', GottchaGenusParser(processor))
        m.add('karp', KarpParser(processor))
        m.add('kslam', KslamParser(processor))
        m.add('megablast', MegablastParser(processor))
        m.add('metaothello', MetaothelloParser(processor))
        m.add('metaphlan2', MetaphlanParser(processor))
        m.add('mmseqs2', Mmseqs2Parser(processor))
        m.add('motus2', Motus2Parser(processor))
        m.add('pathseq', PathseqParser(processor))
        m.add('prophyle', ProphyleParser(processor))
        m.add('taxmaps', TaxmapsParser(processor))
        return m


class ReportParser:
    MATCH_GLOB = None
    IGNORED_NAMES = []
    def __init__(self, processor):
        self.processor = processor
        self.db = processor.db

    def match_filepath(self, filepath):
        if self.MATCH_RE:
            mo = self.MATCH_RE.match(filepath.name)
            return bool(mo), mo
        return False, ''

    def parse(self, in_f, dbase):
        raise NotImplementedError

    def fn_prefix(self, filepath):
        basename = filepath.name

        suff_len = len(self.MATCH_GLOB[1:])
        fn_prefix = basename[:-1*suff_len]
        return fn_prefix


class MetaphlanParser(ReportParser):
    MATCH_GLOB = '*.metaphlan2.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.metaphlan2\.tsv')

    def parse(self, in_f, dbase):
        try:
            next(in_f)
        except StopIteration:
            return
        ranks = {
            'k': 'kingdom',
            'p': 'phylum',
            'c': 'class',
            'o': 'order',
            'f': 'family',
            'g': 'genus',
            's': 'species',
            't': 'subspecies',
            }
        for row in csv.reader(in_f, delimiter='\t'):
            d = dbase.copy()
            taxid = None
            if row[0] == 'unclassified':
                name = 'unclassified'
            else:
                name_list = row[0].split('|')
                deepest = name_list[-1]
                rank_char, name = deepest.split('__')
                d['rank'] = ranks[rank_char]

                res = self.find_name(name)
                if res:
                    name, taxid = res

            abundance = float(row[1]) / 100
            classified_count = self.processor.classified_counts[(dbase['sample'], 'metaphlan2', 'v20')]
            total_count = self.processor.total_reads[(dbase['sample'])]
            classified_prop = classified_count / total_count
            d.update({
                'classifier': 'metaphlan2',
                'database': 'v20',
                'cum_abundance': float(row[1]) / 100 * classified_prop,
                'name': name,
                'taxid': taxid,
            })
            if name.startswith('GCF') or name.startswith('PRJNA') or name.endswith('unclassified') or name.startswith('GCA') or name.endswith('noname'):
                d['no_taxid'] = True

            yield d

        yield self.processor.create_unclassified_taxon(d)

    def names(self, name):
        names = []
        names.append(name)

        name = re.sub('_', ' ', name)
        name = re.sub(' sp ', ' sp. ', name)

        names.append(name)

        name = re.sub(' ([IVX]+) ', r' \1. ', name)
        names.append(name)
        return names


    def find_name(self, name):
        for try_name in self.names(name):
            taxid = self.processor.name_to_taxid.get(try_name)
            if taxid is not None:
                return try_name, taxid


class KrakenParser(ReportParser):
    MATCH_GLOB = '*.kraken.*.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.kraken(\.(?P<database>[^.]+))?\.txt')
    CLASSIFIER = 'kraken'
    is_paired_count = True

    def __init__(self, *args, **kwargs):
        self.fix_species_group = kwargs.pop('fix_species_group', False)
        super().__init__(*args, **kwargs)

    def match_filepath(self, filepath):
        if self.MATCH_RE:
            mo = self.MATCH_RE.match(filepath.name)
            if not mo:
                return False, ''
            elif 'bracken' in mo.groupdict().get('database', ''):
                return False, 'do_not_parse'
            return True, mo
        return False, ''

    def parse(self, in_f, dbase, classifier=None, database='default'):
        classifier = classifier or self.CLASSIFIER

        total_reads = self.processor.total_reads[dbase['sample']]
        for row in csv.reader(in_f, delimiter='\t'):
            taxid = int(row[4])
            rank = None
            # Artifact of Kraken labelling species group as species
            if row[3] == 'S' and self.fix_species_group:
                rank = self.db.ranks.get(taxid)
                if rank and rank == 'species group':
                    rank = 'species group'
            rank = rank or KRAKEN_CODE_BACK.get(row[3])
            cum_reads = int(row[1])
            if dbase['paired'] and self.is_paired_count:
                sample = dbase['sample']
                if not (classifier == 'diamond' and sample.endswith('270')):
                    cum_reads *= 2
            d = dbase.copy()
            d.update({
                'classifier': classifier,
                'database': database,
                'cum_reads': cum_reads ,
                'cum_abundance': cum_reads / total_reads,
                'rank': rank,
                'name': row[5].lstrip(),
                'taxid': taxid,
            })
            yield d


class Kraken2Parser(KrakenParser):
    MATCH_GLOB = '*.kraken2.*.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.kraken2(\.(?P<database>[^.]+))?\.txt')
    CLASSIFIER = 'kraken2'


class BrackenParser(ReportParser):
    MATCH_GLOB = '*.bracken.*.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.bracken(\.(?P<database>[^.]+))?\.txt')
    CLASSRANK = 'species'

    def parse(self, in_f, dbase, classifier='bracken', database='default', fix_species_group=False):
        next(in_f)

        classified_reads = 0
        for row in csv.reader(in_f, delimiter='\t'):
            name = row[0]
            taxid = int(row[1])
            rank_code = row[2]
            rank = KRAKEN_CODE_BACK.get(rank_code)
            kraken_reads = int(row[3])
            new_reads = int(row[5])
            if dbase['paired']:
                new_reads *= 2

            classified_reads += new_reads
            fraction = float(row[6])

            # rank = None
            # if rank_code == 'S' and fix_species_group:
            #     db_rank = self.db.ranks.get(taxid)
            #     if rank and rank == 'species group':
            #         rank = 'species group'

            d = dbase.copy()
            total_reads = self.processor.total_reads[d['sample']]
            d.update({
                'classifier': classifier,
                'database': database,
                'cum_reads': new_reads,
                'cum_abundance': new_reads / total_reads,
                'rank': rank,
                'name': name,
                'taxid': taxid,
            })
            yield d
        yield {
            'sample': d['sample'],
            'classifier': classifier,
            'database': database,
            'cum_reads': total_reads - classified_reads,
            'cum_abundance': (total_reads - classified_reads) / total_reads,
            'name': 'unclassified',
            'taxid': 0,
        }



class BrackenGenusParser(BrackenParser):
    MATCH_GLOB = '*.bracken_genus.*.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.bracken_genus(\.(?P<database>[^.]+))?\.txt')
    CLASSRANK = 'genus'


class KrakenhllParser(ReportParser):
    MATCH_GLOB = '*.krakenhll.*.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.krakenhll\.(?P<database>[^.]+).txt')

    def parse(self, in_f, dbase, database='default', fix_species_group=False):
        total_reads = self.processor.total_reads[dbase['sample']]
        next(in_f)
        ds = []
        kmers_sum = collections.Counter()

        for row in csv.reader(in_f, delimiter='\t'):
            taxid = int(row[6])
            # rank = None
            # if row[7] == 'species' and fix_species_group:
            #     rank = self.db.ranks.get(taxid)
            #     if rank and rank == 'species group':
            #         rank = 'species group'
            # rank = rank or KRAKEN_CODE_BACK.get(row[7])
            rank = row[7]
            kmers = int(row[3])
            uniq_reads = int(row[2])
            cum_reads = int(row[1])
            if dbase['paired']:
                cum_reads *= 2
                uniq_reads *= 2

            d = dbase.copy()
            d.update({
                'classifier': 'krakenhll',
                'database': database,
                'cum_reads': cum_reads,
                'unique_reads': uniq_reads,
                'kmers': kmers,
                'cum_abundance': cum_reads / total_reads,
                'rank': rank,
                'name': row[8].lstrip(),
                'taxid': taxid,
            })
            if rank == 'species':
                kmers_sum['species'] += kmers
            elif rank == 'genus':
                kmers_sum['genus'] += kmers
            ds.append(d)
        for d in ds:
            if d['rank'] == 'species':
                d['cum_abundance'] = d['cum_abundance']
            elif d['rank'] == 'genus':
                d['cum_abundance'] = d['cum_abundance']
            yield d


class ProphyleParser(KrakenParser):
    MATCH_GLOB = '*.prophyle.*.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.prophyle(\.(?P<database>[^\.]*))?\.txt')
    CLASSIFIER = 'prophyle'


class MegablastParser(KrakenParser):
    MATCH_GLOB = '*.megablast.*.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.megablast\.(?P<database>[^.]*)\.txt')
    CLASSIFIER = 'megablast'
    is_paired_count = False


class DiamondKrakenParser(KrakenParser):
    MATCH_GLOB = '*.diamondkraken.report'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.diamondkraken(\.(?P<database>[^\.]*))?\.report')
    CLASSIFIER = 'diamond'
    is_paired_count = False


class KaijuParser(ReportParser):
    MATCH_GLOB = '*.kaiju.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.kaiju\.(?P<database>[^.]*)\.tsv')
    CLASSRANK = 'species'
    IGNORED_NAMES = ['classified_over_species']

    def parse(self, in_f, dbase, database=None, abundance_col='species_abundance'):
        try:
            next(in_f)
        except StopIteration:
            return
        for row in csv.reader(in_f, delimiter='\t'):
            d = dbase.copy()
            taxid = row[3]
            if taxid == 'NA':
                taxid = 0
            else:
                taxid = int(taxid)
            abundance = float(row[1]) / 100
            d.update({
                'classifier': 'kaiju',
                'database': database,
                'cum_abundance': abundance,
                'taxid': taxid,
                'full_name': row[4],
            })
            d[abundance_col] = abundance
            full_name_parts = d['full_name'].split(';')

            if len(full_name_parts) > 1:
                d['name'] = full_name_parts[-2].strip()
            else:
                d['name'] = full_name_parts[0].strip()

            if d['name'] == 'cannot be assigned to a species':
                del d[abundance_col]
                d['name'] = 'classified_over_species'
            yield d

class KaijuGenusParser(KaijuParser):
    MATCH_GLOB = '*.kaiju_genus.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.kaiju_genus\.(?P<database>[^.]*)\.tsv')
    CLASSRANK = 'genus'
    IGNORED_NAMES = ['classified_over_genus']

    def parse(self, in_f, dbase, database=None, abundance_col='genus_abundance'):
        for d in super().parse(in_f, dbase, database=database, abundance_col=abundance_col):
            if d['name'] == 'cannot be assigned to a genus':
                d.pop(abundance_col)
                d['name'] = 'classified_over_genus'
                d['name'] = 'unclassified'
            yield d

class TaxmapsParser(ReportParser):
    MATCH_GLOB = '*.taxmaps.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.taxmaps(\.(?P<database>[^\.]*))?\.tsv')
    IGNORED_NAMES = ['Filtered out']

    def parse(self, in_f, dbase, database=None):
        species_abund_sum = 0
        genus_abund_sum = 0
        ds = []
        for row in csv.reader(in_f, delimiter='\t'):
            d = dbase.copy()
            d.update({
                'classifier': 'taxmaps',
                'database': database,
                'rank': row[1],
                'name': row[4],
                'cum_abundance': float(row[6]) / 100,
            })

            if d['rank'] == 'species':
                species_abund_sum += d['cum_abundance']
            # Filtered out taxid
            if row[0] != '-':
                d['taxid'] = int(row[0])
            ds.append(d)
        for d in ds:
            if d['rank'] == 'species':
                d['species_abundance'] = d['cum_abundance'] / species_abund_sum
            yield d

class CentkrakenParser(KrakenParser):
    MATCH_GLOB = '*.centkraken.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.centkraken(\.(?P<database>[^\.]*))?\.tsv')
    CLASSIFIER = 'centrifuge'


class CentrifugeParser(ReportParser):
    MATCH_GLOB = '*.centrifuge.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.centrifuge(\.(?P<database>[^\.]*))?\.tsv')

    def __init__(self, processor, active=False, *args, **kwargs):
        self.active = active
        super().__init__(processor, *args, **kwargs)

    def parse(self, in_f, dbase, database=None):
        try:
            next(in_f)
        except StopIteration:
            return
        for row in csv.reader(in_f, delimiter='\t'):
            taxid = int(row[1])
            if self.db.ranks[taxid] == 'species group':
                rank = 'species group'
            else:
                rank = row[2]

            d = dbase.copy()
            d.update({
                'classifier': 'centrifuge',
                'rank': rank,
                # 'unique_reads': row[5],
                'name': row[0],
                'taxid': int(row[1]),
                'cls_cum_abundance': float(row[6]),
                'species_abundance': float(row[6]),
            })
            yield d

class GottchaParser(ReportParser):
    MATCH_GLOB = '*.gottcha.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.gottcha\.tsv')
    CLASSRANK = 'species'

    def parse(self, in_f, dbase):
        try:
            next(in_f)
        except StopIteration:
            return
        rows = collections.defaultdict(float)
        for row in csv.reader(in_f, delimiter='\t'):
            d = dbase.copy()
            d.update({
                'classifier': 'gottcha',
                'rank': row[0],
                'name': row[1],
                'bp': int(row[4]),
            })
            rows[d['name']] = d
        total_bps = sum(x['bp'] for x in rows.values())
        for name, d in rows.items():
            d['cls_cum_abundance'] = d['bp'] / total_bps
            d['cum_abundance'] = d['bp'] / total_bps
            if name.lower().startswith('unassigned'):
                d['no_taxid'] = True
            yield d

class GottchaGenusParser(GottchaParser):
    MATCH_GLOB = '*.gottcha_genus.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.gottcha_genus\.tsv')
    CLASSRANK = 'genus'

class KarpParser(ReportParser):
    MATCH_GLOB = '*.karp.nofilter.collapse.freqs'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.karp\.nofilter\.collapse\.freqs')

    def parse(self, in_f, dbase):
        try:
            next(in_f)
        except StopIteration:
            return
        rows = collections.defaultdict(float)

        sample = dbase['sample']
        total_reads = self.processor.total_reads[sample]
        total_e_reads = 0
        for row in csv.reader(in_f, delimiter='\t'):
            abund = float(row[4])
            if row[2] == 'Root':
                taxid = 1
                total_e_reads = abund
            else:
                taxid = int(row[2])
            d = dbase.copy()
            d.update({
                'classifier': 'karp',
                'cum_abundance': abund / total_reads,
                'reads': abund,
                'taxid': taxid,
            })
            yield d
        yield {
            'sample': sample,
            'classifier': 'karp',
            'cum_abundance': (total_reads - total_e_reads) / total_reads,
            'name': 'unclassified',
            'taxid': 0,
        }

class PathseqParser(ReportParser):
    MATCH_GLOB = '*.pathseq.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.pathseq\.txt')

    def parse(self, in_f, dbase):
        try:
            next(in_f)
        except StopIteration:
            return
        rows = {}
        total_reads = self.processor.total_reads[dbase['sample']]
        norm_score_sum = 0
        ds = []
        for row in csv.reader(in_f, delimiter='\t'):
            taxid = int(row[0])
            # Not unique or cumulative. Can count multiple times as multihit
            reads = float(row[7])
            score = float(row[5])
            if taxid == 1:
                root_score = score
            d = dbase.copy()
            d.update({
                'classifier': 'pathseq',
                'rank': row[2],
                'name': row[3],
                'score': score,
                'taxid': taxid,
                'norm_score': reads,
                'reads': reads,
                'database': 'default',
            })
            ds.append(d)
        for d in ds:
            # Let's try without using total_reads
            d['cum_abundance'] = d['score'] / total_reads
            # d['cum_abundance'] = reads / total_reads
            yield d
        u_abund = (total_reads - root_score) / total_reads
        yield self.processor.create_unclassified_taxon(
            d, unclassified_abundance=u_abund)

class MetaothelloParser(ReportParser):
    MATCH_GLOB = '*.metaothello.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.metaothello(\.(?P<database>[^\.]*))?\.tsv')

    def parse(self, in_f, dbase, database=None):
        for row in csv.reader(in_f, delimiter='\t'):
            rank = row[3]
            if rank == 'species':
                classrank = 'species'
            elif rank == 'genus':
                classrank = 'genus'
            else:
                log.error('classrank missing: %s', row)
            if row[1] == 'unclassified':
                row[3] = ''
            d = dbase.copy()
            d.update({
                'classifier': 'metaothello',
                'database': database,
                'rank': row[3],
                'name': row[1],
                'taxid': row[0],
                'classrank': classrank,
                'cum_abundance': float(row[2])
            })
            yield d

class Mmseqs2Parser(ReportParser):
    MATCH_GLOB = '*.mmseqs2.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.mmseqs2(\.(?P<database>[^\.]*))?\.tsv')

    def parse(self, in_f, dbase, database=None):
        for row in csv.reader(in_f, delimiter='\t'):
            rank = row[3]
            if rank == 'species':
                classrank = 'species'
            elif rank == 'genus':
                classrank = 'genus'
            else:
                log.error('classrank missing: %s', row)
            if row[1] == 'unclassified':
                row[3] = ''
            d = dbase.copy()
            d.update({
                'classifier': 'mmseqs2',
                'database': database,
                'rank': row[3],
                'name': row[1],
                'taxid': row[0],
                'classrank': classrank,
                'cum_abundance': float(row[2])
            })
            yield d


class Motus2Parser(ReportParser):
    MATCH_GLOB = '*.motus.txt'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.motus\.txt')

    def parse(self, in_f, dbase, database=None):
        taxa = []
        cum_reads = collections.Counter()
        total_reads = self.processor.total_reads[dbase['sample']]
        cls_reads = 0
        for row in csv.reader(in_f, delimiter='\t'):
            if row[0].startswith('#'):
                continue
            name = row[1]
            name = re.sub(r' \[C\]', '',  name)
            try:
                taxid = int(row[2])
            except:
                continue
            num_reads = float(row[3])
            if num_reads == 0:
                continue

            cum_reads[taxid] += num_reads
            path = self.db.parent_path(taxid, cache=self.processor.parent_path_cache, warn_missing=False)
            for taxid in path:
                cum_reads[taxid] += num_reads
                cls_reads += num_reads

        for taxid, cum_reads in cum_reads.items():
            rank = self.db.ranks.get(taxid)
            names = self.db.names.get(taxid)[0]
            if not names:
                continue
            name = names[0]
            d = {
                'sample': dbase['sample'],
                'classifier': 'motus2',
                'database': 'default',
                'rank': rank,
                'name': name,
                'taxid': taxid,
                'classrank': 'all',
                'cum_abundance': cum_reads / total_reads,
            }
            yield d
        u_abund = (total_reads - cls_reads) / total_reads

        d = dbase.copy()
        d['classifier'] = 'motus2'
        d['database'] = 'default'
        yield self.processor.create_unclassified_taxon(d, unclassified_abundance=u_abund)


class KslamParser(ReportParser):
    MATCH_GLOB = '*.kslam.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.kslam(\.(?P<database>[^\.]*))?\.tsv')

    def parse(self, in_f, dbase, database=None):
        taxa = []

        total_reads = self.processor.total_reads[dbase['sample']]
        classified_count = self.processor.classified_counts[(dbase['sample'], 'kslam', database)]
        classified_prop = classified_count / total_reads
        missing_taxons = []
        for row in csv.reader(in_f, delimiter='\t'):
            abundance_percent = float(row[1])
            if row[0] == '':
                # This is normal if k-SLAM can't find name for a taxid - probably merged
                missing_taxons.append(abundance_percent)

            name = row[0]
            taxid = self.processor.name_lower_to_taxid.get(name.lower())
            name = self.processor.name_lower_to_original.get(name, name)
            d = dbase.copy()
            d.update({
                'classifier': 'kslam',
                'database': database,
                'name': name,
                'taxid': taxid,
                'unique_abundance': abundance_percent / 100 * classified_prop
            })
            if taxid is not None:
                d['taxid'] = int(taxid)
            taxa.append(d)

        if missing_taxons:
            log.info('K-SLAM: missing taxon names for file: %s, abundances %s', in_f.name, missing_taxons)
        cum_abunds = collections.Counter()
        for d in taxa:
            if d['name'] == '':
                continue

            if 'taxid' in d:
                taxid = d['taxid']
                d['rank'] = self.db.ranks.get(taxid)
                cum_abunds[taxid] += d['unique_abundance']
                path = self.db.parent_path(taxid, cache=self.processor.parent_path_cache, warn_missing=False)
                for taxid in path:
                    cum_abunds[taxid] += d['unique_abundance']

        seen_taxa = set()
        for d in taxa:
            if 'taxid' in d:
                seen_taxa.add(d['taxid'])
                d['cum_abundance'] = cum_abunds[d['taxid']]

        for taxid, cum_abund in cum_abunds.items():
            if taxid in seen_taxa:
                continue
            rank = self.db.ranks.get(taxid)
            d = dbase.copy()
            d['rank'] = rank
            d['database'] = database
            d['taxid'] = taxid
            d['name'] = self.db.names[taxid][0]
            d['classifier'] = 'kslam'
            d['cum_abundance'] = cum_abund
            d['classrank'] = 'all'
            taxa.append(d)
        taxa.append(self.processor.create_unclassified_taxon(d, unclassified_abundance=1 - classified_prop + sum(missing_taxons)))
        return taxa


class ClarkParser(ReportParser):
    MATCH_GLOB = '*.clark.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.clark(\.(?P<database>[^.]*))?\.csv')
    CLASSIFIER = 'clark'
    CLASSRANK = 'species'

    def parse(self, in_f, dbase, database=None):
        try:
            next(in_f)
        except StopIteration:
            return
        for row in csv.reader(in_f, delimiter=','):
            taxid = row[1]
            rank = self.CLASSRANK
            if taxid == 'UNKNOWN':
                taxid = 0
                rank = ''
            d = dbase.copy()
            d.update({
                'classifier': self.CLASSIFIER,
                'database': database,
                'rank': rank,
                'name': row[0],
                'taxid': taxid,
                'cum_abundance': float(row[4]) / 100,
            })
            yield d

class ClarkSParser(ClarkParser):
    MATCH_GLOB = '*.clark_s.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.clark_s(\.(?P<database>[^.]*))?\.csv')
    CLASSIFIER = 'clark_s'
    CLASSRANK = 'species'


class ClarkGenusParser(ClarkParser):
    MATCH_GLOB = '*.clark_genus.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.clark_genus(\.(?P<database>[^.]*))?\.csv')
    CLASSIFIER = 'clark'
    CLASSRANK = 'genus'

class ClarkSGenusParser(ClarkSParser):
    MATCH_GLOB = '*.clark_s_genus.*.tsv'
    MATCH_RE = re.compile(r'(?P<fn_prefix>.*)\.clark_s_genus(\.(?P<database>[^.]*))?\.csv')
    CLASSIFIER = 'clark_s'
    CLASSRANK = 'genus'


def read_classified_counts(cc_f):
    classified_counts = {}
    for row in csv.reader(cc_f, delimiter='\t'):
        name, count = row[0], row[1]
        count = int(count)
        if '.refseqc' in name:
            db = 'refseqc'
        elif '.bvp' in name:
            db = 'bvp'
        elif '.nr' in name:
            db = 'nr'
        elif 'metaphlan2' in name:
            db = 'v20'
        elif 'pathseq' in name:
            db = 'default'
        elif 'kslam' in name:
            db = 'default'
        else:
            raise Exception


        if 'metaphlan2' in name:
            sample = name[:name.find('metaphlan2') - 1]
            classifier = 'metaphlan2'
        elif 'pathseq' in name:
            sample = name[:name.find('pathseq') - 1]
            classifier = 'pathseq'
        elif 'kslam' in name:
            sample = name[:name.find('kslam') - 1]
            classifier = 'kslam'


        if sample.endswith('.art'):
            sample = sample[:-4]

        classified_counts[(sample , classifier, db)] = count
    return classified_counts


def compile_reports(args):
    assert args.files or args.dir, "One of --files or --dir must be specified"
    conf = load_config(args.config)

    db = ncbitax.TaxonomyDb(tax_dir=args.tax_dir, load_nodes=True, load_names=True, load_merged=True,
                            scientific_names_only=False)

    with open(args.output, 'wt') as out_f:

        writer = ReportTableWriter(out_f)
        proc = ReportProcessor(db, conf)

        proc.total_reads = {}
        if args.read_counts:
            with open(args.read_counts) as reads_f:
                for row in csv.reader(reads_f, delimiter='\t'):
                    sample, count = row
                    if sample.endswith('.art'):
                        sample = sample[:-4]
                    proc.total_reads[sample] = int(count)
        if args.classified_counts:
            with open(args.classified_counts) as cc_f:
                proc.classified_counts = read_classified_counts(cc_f)

        taxa = []
        if args.files:
            for fn in args.files:
                basename = os.path.basename(fn)
                res = proc.process_file(basename, open(fn))
                if res:
                    taxa.extend(res)
                #sys.exit()
        elif args.dir:
            for x in Path(args.dir).iterdir():
                if x.is_dir():
                    continue
                basename = x.name
                with x.open() as in_f:
                    res = proc.process_file(basename, in_f)
                    if res:
                        taxa.extend(res)
        proc.sum_rank_abundances(taxa)
        for d in taxa:
            writer.write(d)
        with open(args.rank_abundances, 'wt') as out_f:
            for k, v in proc.rank_abundances.items():
                parts = [str(x) for x in list(k) + [v]]
                print('\t'.join(parts), file=out_f)

        if args.missing_parents:
            with open(args.missing_parents, 'wt') as f:
                missing_parents = {k for k, v in proc.parent_path_cache.items() if v is None}
                log.info('Missing parents for %s taxids', len(missing_parents))
                for taxid in sorted(list(missing_parents)):
                    if taxid == 0:
                        continue
                    print(taxid, file=f)



def add_command(subparsers):
    parser = subparsers.add_parser('compile', description='Compile a bunch of report files together')
    parser.add_argument('output')
    parser.add_argument('--config')
    parser.add_argument('--tax-dir', required=True)
    parser.add_argument('--classified-counts', required=True)
    parser.add_argument('--read-counts', required=True)
    parser.add_argument('--rank-abundances', required=True)
    parser.add_argument('--missing-parents')
    parser.add_argument('--files', nargs='+')
    parser.add_argument('--dir')
    parser.set_defaults(func=compile_reports)
