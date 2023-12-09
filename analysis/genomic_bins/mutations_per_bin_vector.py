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
    sum_probs_dict = defaultdict(lambda: defaultdict(float))
    signatures = list()
    with open(muts_f, 'r') as fd:
        for i, line in enumerate(fd):
            # Parse data
            if i > 0:
                data = line.strip().split('\t')
                chrom, pos = data[1:3]
                probabilities = data[4:]
                for interval in tree[chrom][int(pos)]:
                    binid = interval.data
                    # Load
                    for s, p in list(zip(signatures, probabilities)):
                        sum_probs_dict[binid][s] += float(p)
            # Parse header
            else:
                header = line.strip().split('\t')
                signatures = sorted(header[4:])

    # Save
    dict_to_json = {}
    for binid, sigs_to_muts in sum_probs_dict.items():
        dict_to_json[binid] = {}
        for sig, mutations in sigs_to_muts.items():
            dict_to_json[binid][sig] = mutations
    with open(output_f, 'w') as ofd:
        json.dump(dict_to_json, ofd)


if __name__ == '__main__':
    main()
