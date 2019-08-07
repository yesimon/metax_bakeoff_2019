#!/usr/bin/env python3
import argparse
import itertools
import collections
import logging
import heapq
import os
import sys
import gzip

import cigar
import numpy as np
import pyfaidx
import pysam

# log = logging.getLogger('merge_sam')
log = logging.getLogger()

def guess_sam_encoding(sam_fn, mode='r'):
    '''Use filename extension to guess binary bam vs text sam'''
    if sam_fn.endswith('.sam'):
        pass
    elif sam_fn.endswith('.bam'):
        mode += 'b'
    else:
        print('Unknown file extension for sam file: {}, assuming bam'.format(os.path.splitext(sam_fn), file=sys.stderr))
        mode += 'b'
    return mode


def alignment_score(sam1, mate_sam1=None):
    if sam1.has_tag('AS'):
        read_as = sam1.get_tag('AS')
    else:
        read_as = 0

    if mate_sam1:
        if mate_sam1.has_tag('AS'):
            mate_as = mate_sam1.get_tag('AS')
            read_as += mate_as
    return read_as


def fullmerge(*iterables, key=None):
    giters = []
    for it in iterables:
        giter = itertools.groupby(it, key)
        giters.append(giter)

    while True:
        try:
            glist = []
            for pair in giters:
                gkey, giter = next(pair)
                glist.extend(list(giter))
            yield glist
        except StopIteration:
            return


def fullmerge_parts(*iterables, key=None):
    giters = []
    for it in iterables:
        giter = itertools.groupby(it, key)
        giters.append(giter)

    while True:
        try:
            glist = []
            for i, pair in enumerate(giters):
                gkey, giter = next(pair)
                glist.append(list(giter))
            yield glist
        except StopIteration:
            return

def imerge(*iterables):
    '''Merge multiple sorted inputs into a single sorted output.

    Equivalent to:  sorted(itertools.chain(*iterables))

    >>> list(imerge([1,3,5,7], [0,2,4,8], [5,10,15,20], [], [25]))
    [0, 1, 2, 3, 4, 5, 5, 7, 8, 10, 15, 20, 25]

    '''
    heappop, siftup, _StopIteration = heapq.heappop, heapq._siftup, StopIteration

    h = []
    h_append = h.append
    for it in map(iter, iterables):
        try:
            next = it.next
            h_append([next(), next])
        except _StopIteration:
            pass
    heapq.heapify(h)

    while 1:
        try:
            while 1:
                v, next = s = h[0]      # raises IndexError when h is empty
                yield v
                s[0] = next()           # raises StopIteration when exhausted
                siftup(h, 0)            # restore heap condition
        except _StopIteration:
            heappop(h)                  # remove empty iterator
        except IndexError:
            return

def read_pairs(sam1s):
    for i in range(0, len(sam1s), 2):
        yield (sam1s[i], sam1s[i+1])

def pair_alns_as(sam1_groups):
    best_a_score = 0
    best_sams = None
    for sam1_group in sam1_groups:
        primary_alns = [sam1 for sam1 in sam1_group if not sam1.flag & 0x900 and sam1.flag & 0xc != 0xc]
        if not primary_alns:
            continue
        if len(primary_alns) == 1:
            a_score = alignment_score(primary_alns[0])
            if a_score > best_a_score:
                best_a_score = a_score
                best_sams = primary_alns
        elif len(primary_alns) == 2:
            a_score = alignment_score(primary_alns[0], primary_alns[1])
            if a_score > best_a_score:
                best_a_score = a_score
                best_sams = primary_alns
        else:
            raise ValueError(sam1_group)
    if not best_sams:
        return
    for sam1_group in sam1_groups:
        for sam1 in sam1_group:
            if sam1 in best_sams or sam1.flag & 0xc == 0xc:
                continue
            sam1.flag |= 0x100

def add_command(subparsers):
    parser = subparsers.add_parser('merge-sam')
    parser.add_argument('sams', nargs='+')
    parser.add_argument('-o', '--output')
    parser.add_argument('-l', '--loglevel', default=logging.WARNING)
    parser.set_defaults(func=merge_sam)


def merge_sam(args):
    if args.loglevel:
        loglevel = getattr(logging, args.loglevel.upper())
        print(args.loglevel)
        print(loglevel)
        assert type(loglevel) is int
        log.setLevel(loglevel)

    in_sams = []
    for sam in args.sams:
        if os.path.isfile and os.stat(sam).st_size <= 71:
            continue
        in_sams.append(pysam.AlignmentFile(sam, guess_sam_encoding(sam, 'r')))


    header = in_sams[0].header
    all_sq = [sq['SN'] for sq in header['SQ']]
    tid_translation = {sq['SN']: i for i, sq in enumerate(all_sq)}
    all_sq = set(all_sq)
    sq_len = len(all_sq)
    pg_ids = collections.Counter(pg['ID'] for pg in header['PG'])
    log.info('Read first headers')
    for sam in in_sams[1:]:
        new_header_sq = []
        for i, sq in enumerate(sam.header['SQ']):
            sn = sq['SN']

            if sn not in all_sq:
                tid_translation[sn] = sq_len
                new_header_sq.append(sq)
                all_sq.add(sn)
                sq_len += 1
        for pg in sam.header['PG']:
            pid = pg['ID']
            if pid in pg_ids:
                pg['ID'] = '{}_{}'.format(pid, pg_ids[pid])
                pg_ids[pid] += 1
            header['PG'].append(pg)

        header['SQ'].extend(new_header_sq)

    log.info('Reading headers finished')
    out_sam = pysam.AlignmentFile(args.output, guess_sam_encoding(args.output, 'w'), header=header)
    for i, sam1_groups in enumerate(fullmerge_parts(*in_sams, key=lambda x: x.query_name)):
        pair_alns_as(sam1_groups)
        sam1_group = list(itertools.chain(*sam1_groups))

        first_aligned = any(sam1 for sam1 in sam1_group if sam1.flag & 0x40 and not sam1.flag & 0x4)
        second_aligned = any(sam1 for sam1 in sam1_group if sam1.flag & 0x80 and not sam1.flag & 0x4)

        first_seen = False
        second_seen = False
        for sam1 in sam1_group:
            if sam1.flag & 0xc == 0xc:
                if first_aligned or second_aligned:
                    continue
                if sam1.flag & 0x40:
                    if first_seen:
                        continue
                    else:
                        first_seen = True
                elif sam1.flag & 0x80:
                    if second_seen:
                        continue
                    else:
                        second_seen = True

            # if sam1.flag & 0x40 and sam1.flag & 0x4:
            #     if first_aligned or first_seen:
            #         continue
            #     else:
            #         if first_aligned:
            #             sam1.flag -= 0x4
            #         if second_aligned and sam1.flag & 0x8:
            #             sam1.flag -= 0x8
            #         first_seen = True
            # elif sam1.flag & 0x80 and sam1.flag & 0x4:
            #     if second_aligned or second_seen:
            #         continue
            #     else:
            #         if second_aligned:
            #             sam1.flag -= 0x4
            #         if first_aligned and sam1.flag & 0x8:
            #             sam1.flag -= 0x8
            #         second_seen = True
            if sam1.reference_id >= 0:
                new_tid = tid_translation.get(sam1.reference_name)
                if new_tid:
                    sam1.reference_id = new_tid
            if sam1.next_reference_id >= 0:
                new_tid = tid_translation.get(sam1.next_reference_name)
                if new_tid:
                    sam1.next_reference_id = new_tid
            out_sam.write(sam1)
