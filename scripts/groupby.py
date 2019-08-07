#!/usr/bin/env python3
import argparse
import sys
import collections


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='/dev/stdout')
    parser.add_argument('--input', default='/dev/stdin')
    parser.add_argument('--action', default='count')
    parser.add_argument('--stats')

    args = parser.parse_args()

    totals = collections.Counter()
    lengths = collections.defaultdict(int)
    for line in sys.stdin:
        line = line.strip()
        parts = line.split('\t')
        name = parts[0]
        cov = int(parts[2])
        lengths[name] = int(parts[1])
        totals[name] += cov
    for name, count in totals.items():
        avg = count / lengths[name]
        print('{}\t{}'.format(name, avg))




if __name__ == '__main__':
    main()
