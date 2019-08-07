#!/usr/bin/env python3
import argparse
import gzip
import bz2
from Bio import SeqIO


def compressed_open(fname, *args, **kwargs):
    if fname.endswith('.gz'):
        return gzip.open(fname, *args, **kwargs)
    elif fname.endswith('.bz2'):
        return bz2.open(fname, *args, **kwargs)
    else:
        return open(fname, *args, **kwargs)



def main():
    parser = argparse.ArgumentParser(description='Fix problems with raw fastqs')
    parser.add_argument('inputs', nargs='+', help='Input paired fastq files.')
    parser.add_argument('-o', '--output', nargs='+', help='Output paired fastq files.')
    parser.add_argument('--unpaired-output', nargs='+', help='Output paired fastq files.')
    args = parser.parse_args()

    assert len(args.inputs) == 2
    f_ids = [set(), set()]
    for i in range(2):
        for rec in SeqIO.parse(compressed_open(args.inputs[i], 'rt'), 'fastq'):
            if len(rec.seq) == 0:
                continue
            if rec.id.endswith('/1') or rec.id.endswith('/2'):
                f_ids[i].add(rec.id[:-2])

    both_ids = f_ids[0] & f_ids[1]

    for i in range(2):
        with compressed_open(args.output[i], 'wt') as f:
            for rec in SeqIO.parse(compressed_open(args.inputs[i], 'rt'), 'fastq'):
                seq_id = rec.id
                if rec.id.endswith('/1') or rec.id.endswith('/2'):
                    seq_id = rec.id[:-2]
                if seq_id in both_ids:
                    SeqIO.write(rec, f, 'fastq')


if __name__ == '__main__':
    main()
