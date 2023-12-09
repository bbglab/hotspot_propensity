from collections import defaultdict
import json

import click
from intervaltree import IntervalTree
import pandas as pd

autosomes = list(map(lambda x: f'chr{x}', range(1, 23)))


@click.command()
@click.option('-h', '--hotspots_f', default=None, required=True)
@click.option('-b', '--bins_f', default=None, required=True)
@click.option('-o', '--output_f', default=None, required=True)
def main(hotspots_f, bins_f, output_f):

    # Load bins in autosomes
    bins_df = pd.read_csv(bins_f, sep='\t', header=0)

    # Load bins into intervaltree
    tree = defaultdict(IntervalTree)
    for binid in bins_df['BINID'].tolist():
        chrom, start_end = binid.split(':')
        chrom = chrom[3:]
        start, end = start_end.split('-')
        tree[chrom].addi(int(start), int(end) + 1, binid)  # +1 open end

    # Read hotspots and intersect
    hotspots_per_bin = defaultdict(lambda: defaultdict(int))
    hotspot_sigs = defaultdict(lambda: defaultdict(float))
    signatures = set()
    with open(hotspots_f, 'r') as fd:
        next(fd)
        for line in fd:
            hotspot, chrom, pos, _, _, signature, _, prob = line.strip().split()
            hotspot_sigs[hotspot][signature] = prob
            signatures.add(signature)
    # Assign hotspot to signature by maximum likelihood
    for hotspot, data in hotspot_sigs.items():
        sorted_sigs = sorted([(s, float(p)) for s, p in data.items()], key=lambda x: x[1], reverse=True)
        best_sig = sorted_sigs[0][0]
        chrom, pos, _ = hotspot.split('_')
        # Intersect with bins
        for interval in tree[str(chrom)][int(pos)]:
            hotspots_per_bin[interval.data][best_sig] += 1

    # Save
    dict_to_json = {}
    for binid, sigs_to_hotspots in hotspots_per_bin.items():
        dict_to_json[binid] = {}
        for sig, hotspots in sigs_to_hotspots.items():
            dict_to_json[binid][sig] = hotspots
    with open(output_f, 'w') as ofd:
        json.dump(dict_to_json, ofd)


if __name__ == '__main__':
    main()
