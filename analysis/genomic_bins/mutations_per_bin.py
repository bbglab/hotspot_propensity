from collections import defaultdict
import json

import click
from intervaltree import IntervalTree
import pandas as pd


@click.command()
@click.option('-m', '--muts_f', default=None, required=True)
@click.option('-b', '--bins_f', default=None, required=True)
@click.option('-o', '--output_f', default=None, required=True)
def main(muts_f, bins_f, output_f):

    # Load bins
    bins_df = pd.read_csv(bins_f, sep='\t', header=0)
    selected_bins = set(bins_df['BINID'].unique())

    # Load bins into intervaltree
    tree = defaultdict(IntervalTree)
    for binid in selected_bins:
        chrom, start_end = binid.split(':')
        start, end = start_end.split('-')
        chrom = chrom[3:]
        tree[chrom].addi(int(start), int(end) + 1, binid)

    # Read mutations and intersect
    mutations_per_bin = defaultdict(lambda: defaultdict(int))
    with open(muts_f, 'r') as fd:
        next(fd)
        for line in fd:
            # Parse data
            _, chrom, pos, _, _, _, _, s, _ = line.strip().split('\t')
            for interval in tree[chrom][int(pos)]:
                binid = interval.data
                mutations_per_bin[binid][s] += 1

    # Save
    dict_to_json = {}
    for binid, sigs_to_muts in mutations_per_bin.items():
        dict_to_json[binid] = {}
        for sig, mutations in sigs_to_muts.items():
            dict_to_json[binid][sig] = mutations
    with open(output_f, 'w') as ofd:
        json.dump(dict_to_json, ofd)


if __name__ == '__main__':
    main()
