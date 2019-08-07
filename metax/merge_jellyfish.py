import argparse
import collections
import subprocess
import textwrap
from pathlib import Path


def merge_jellyfish_dbs(args):
    merging_dbs = collections.defaultdict(list)
    for path in Path(args.dirpath).iterdir():
        if path.match('*.jdb_*'):
            db_name = path.stem
            merging_dbs[db_name].append(path)
    # Look for jellyfish dbs that are split in multiple parts
    jf_names = []
    for db_name, jdbs in merging_dbs.items():
        fna_name = db_name
        jf_name = Path(fna_name + '.jf')
        jf_names.append(jf_name)
        if jf_name.exists():
            continue
        if jf_name.is_symlink():
            jf_name.unlink()
        if len(jdbs) < 2:
            jf_name.symlink_to(fna_name + '.jdb_0')
            continue
        print(db_name, jdbs)
        subprocess.run('jellyfish merge -o {} {}'.format(jf_name, ' '.join([str(x) for x in jdbs])), shell=True)
    for jf_name in jf_names:
        jtsv = jf_name.with_suffix('.jtsv')
        jdone = jf_name.with_suffix('.jdone')
        if jdone.exists():
            continue
        subprocess.run('jellyfish dump -t -c -o {} {}'.format(jtsv, jf_name), shell=True)
        jdone.touch()

description = textwrap.dedent('''\
Merge all jellyfish1 taxon dbs that have been split due to exceeding default max size.
''')


def add_command(subparsers):
    parser = subparsers.add_parser('merge-jellyfish-dbs', description=description)
    parser.add_argument('dirpath', help='Directory of .jdb* databases')
    parser.set_defaults(func=merge_jellyfish_dbs)
