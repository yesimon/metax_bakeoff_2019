#!/usr/bin/env python3
import collections
import os
import subprocess
import shlex
import sys
import tempfile
import argparse
from Bio import SeqIO


def proportional(nseats, votes_dict):
    """assign n seats proportionaly to votes using Hagenbach-Bischoff quota
    :param nseats: int number of seats to assign
    :param votes: iterable of int or float weighting each party
    :result: list of ints seats allocated to each party
    """
    items = votes_dict.items()
    names = [x[0] for x in items]
    votes = [x[1] for x in items]
    quota = sum(votes)/(1.+nseats) #force float
    frac = [vote//quota for vote in votes]
    res = [int(f) for f in frac]
    n = nseats - sum(res) #number of seats remaining to allocate
    if n==0:
        return dict(zip(names, res)) #done
    if n<0:
        return dict(zip(names, [min(x,nseats) for x in res])) # see siamii's comment
    #give the remaining seats to the n parties with the largest remainder
    remainders=[ai-bi for ai,bi in zip(frac,res)]
    limit = sorted(remainders,reverse=True)[n-1]
    #n parties with remainter larger than limit get an extra seat
    for i,r in enumerate(remainders):
        if r>=limit:
            res[i]+=1
            n-=1 # attempt to handle perfect equality
            if n==0:
                return dict(zip(names, res)) #done
    raise #should never happen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('reference')
    parser.add_argument('abundances')
    parser.add_argument('--mode', default='files')
    parser.add_argument('--total-reads', default=2000000)
    parser.add_argument('--read-length', default=150)
    parser.add_argument('--wgsim-random-seed', default=1)
    parser.add_argument('--mutation-rate', type=float)
    parser.add_argument('--output', nargs=2)
    parser.add_argument('--spike-abundances', nargs='*')
    parser.add_argument('--output-seq-abundances', type=argparse.FileType('w'))
    parser.add_argument('--output-seq-fasta', type=argparse.FileType('w'))

    parser.add_argument('--simulator', default='wgsim')
    args = parser.parse_args()

    with open(args.abundances) as f:
        if args.abundances.endswith('.yaml'):
            import yaml
            abundances = yaml.load(f)
        elif args.abundances.endswith('.tsv'):
            import csv
            abundances = {k: float(v) for k, v in csv.reader(f)}

    assert abs(sum(abundances.values()) - 1.0) < 0.000001

    if args.spike_abundances:
        spike_abundances = {}
        for i in range(0, len(args.spike_abundances), 2):
            spike_abundances[args.spike_abundances[i]] = float(args.spike_abundances[i + 1])
        remaining = 1 - sum(spike_abundances.values())
        for k, abund in abundances.items():
            abundances[k] = abund * remaining
        abundances.update(spike_abundances)

    assert abs(sum(abundances.values()) - 1.0) < 0.000001

    total_reads = args.total_reads

    if args.output_seq_abundances:
        seq_abundances_proportions(args.reference, abundances, output_file=args.output_seq_abundances,
                                   output_fasta=args.output_seq_fasta)
        sys.exit()

    for i in args.output:
        with open(i, 'w'):
            pass

    if args.simulator == 'wgsim':
        if args.mode == 'files':
            for fa_fn, prop in abundances.items():
                file_reads = int(total_reads * prop)
                reference_fasta = os.path.join(args.reference, fa_fn)
                prefix1 = '.'.join([args.output[0], fa_fn])
                prefix2 = '.'.join([args.output[1], fa_fn])
                os.mkfifo(prefix1)
                os.mkfifo(prefix2)
                if args.mutation_rate:
                    mut_rate = float(args.mutation_rate)
                else:
                    mut_rate = ''
                try:
                    cmd = 'wgsim -N {num_reads} -S 1 -d 400 -1 {read_length} -2 {read_length} -r {mut_rate} {reference_fasta} {prefix1} {prefix2}'.format(
                        num_reads=file_reads,
                        read_length=args.read_length,
                        mut_rate=mut_rate,
                        prefix1=prefix1,
                        prefix2=prefix2,
                        reference_fasta=reference_fasta)
                    # subprocess.Popen(shlex.split(cmd))
                    subprocess.Popen(cmd, shell=True)
                    subprocess.Popen('cat {prefix1} | tee >(wc -l 1>&2 ) >> {output}'.format(prefix1=prefix1, output=args.output[0]), shell=True, executable='/bin/bash')
                    subprocess.run('cat {prefix2} >> {output}'.format(prefix2=prefix2, output=args.output[1]), shell=True)
                finally:
                    os.remove(prefix1)
                    os.remove(prefix2)
    elif args.simulator == 'art':
        if args.mode == 'files':
            for fa_fn, prop in abundances.items():
                file_nreads = int(total_reads * prop)
                reference_fasta = os.path.join(args.reference, fa_fn)
                prefix = '.'.join([args.output[0], fa_fn])

                seqlens = collections.Counter()
                fas = SeqIO.index(reference_fasta, 'fasta')
                nbases = 0
                for name, seqr in fas.items():
                    seqlen = len(seqr)
                    if seqlen < 1000:
                        continue
                    seqlens[name] = seqlen
                    nbases += seqlen
                seq_nreads = proportional(file_nreads, seqlens)
                for name, seqr in fas.items():
                    if name not in seq_nreads:
                        continue
                    nreads = seq_nreads[name]
                    with tempfile.NamedTemporaryFile('wt') as f:
                        SeqIO.write(seqr, f.file, 'fasta')
                        f.flush()
                        cmd = '~/idi/src/art_src_MountRainier_Linux/art_illumina --noALN --rndSeed 1 -nf 0 -ss HS25 -na -p -sam -na -i {reference_fasta} -l {read_length} -c {num_reads} -m 400 -s 50 -o {prefix}'.format(
                            reference_fasta=f.name,
                            read_length=args.read_length,
                            num_reads=nreads,
                            prefix=prefix)
                        subprocess.run(cmd, shell=True)
                        records = SeqIO.parse('{prefix}1.fq'.format(prefix=prefix), 'fastq')
                        nrecords = len(list(records))
                        if nrecords != nreads:
                            print('Sequence {} - len {} reads mismatch {} out of {} expected'.format(name, seqlens[name], nrecords, nreads), file=sys.stderr)
                        subprocess.run('cat {prefix}1.fq >> {output}'.format(prefix=prefix, output=args.output[0]), shell=True)
                        subprocess.run('cat {prefix}2.fq >> {output}'.format(prefix=prefix, output=args.output[1]), shell=True)
                        os.remove('{prefix}1.fq'.format(prefix=prefix))
                        os.remove('{prefix}2.fq'.format(prefix=prefix))
                        # os.remove('{prefix}.sam'.format(prefix=prefix))


def seq_abundances_counts(total_reads=None):
    file_nreads = int(total_reads * prop)
    reference_fasta = os.path.join(args.reference, fa_fn)
    prefix = '.'.join([args.output[0], fa_fn])

    seqlens = collections.Counter()
    fas = SeqIO.index(reference_fasta, 'fasta')
    nbases = 0
    for name, seqr in fas.items():
        seqlen = len(seqr)
        if seqlen < 1000:
            continue
        seqlens[name] = seqlen
        nbases += seqlen
    seq_nreads = proportional(file_nreads, seqlens)



def seq_abundances_proportions(reference_dir, abundances, output_file=None, output_fasta=None):
    output_file = output_file or sys.stdout
    for fa_fn, abundance in abundances.items():
        reference_fasta = os.path.join(reference_dir, fa_fn)

        seqlens = collections.defaultdict(float)
        # fas = SeqIO.index(reference_fasta, 'fasta')
        nbases = 0
        for seqr in SeqIO.parse(reference_fasta, 'fasta'):
        # for name, seqr in fas.items():

            seqlen = len(seqr)
            if seqlen < 1000:
                continue
            if output_fasta:
                SeqIO.write(seqr, output_fasta, 'fasta')
            seqlens[seqr.description] = seqlen
            nbases += seqlen

        for name, seqr in seqlens.items():
            print('\t'.join(['>{}'.format(name), '{:.3g}'.format(abundance * seqlens[name] / nbases)]), file=output_file)



if __name__ == '__main__':
    main()
