import subprocess
import sys
from pkg_resources import resource_filename


def create_centrifuge_map(args):
    exe = resource_filename(__name__, 'create_centrifuge_map.sh')
    ret = subprocess.run([exe] + sys.argv[2:], check=True)


def add_command(subparsers):
    parser = subparsers.add_parser('create-centrifuge-map')
    parser.add_argument('ncbi_genomes_dir')
    parser.set_defaults(func=create_centrifuge_map)
