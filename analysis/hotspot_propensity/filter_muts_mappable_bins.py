from collections import defaultdict

import click
from intervaltree import IntervalTree
import pandas as pd

autosomes = list(map(lambda n: f'chr{n}', range(1,23)))


@click.command()
@click.option('-i', '--input_f', default=None, help='Input directory')
@click.option('-o', '--output_f', default=None, help='Output file')
@click.option('-b', '--bins_f', default=None, help='Output file')
def main(input_f, output_f, bins_f):
    """Filter mutations in mappable bins"""

    # Load bins and get those in autosomes
    bins_df = pd.read_csv(bins_f, sep='\t', header=0)
    bins_autosom_df = bins_df.loc[bins_df['CHR'].isin(autosomes)].copy()

    # Load bins into intervaltree
    tree = defaultdict(IntervalTree)
    columns = bins_autosom_df.columns
    for _, row in bins_autosom_df.iterrows():
        chrom, start, end, binid = [row[c] for c in columns]
        chrom = chrom[3:]
        tree[chrom].addi(start, end + 1, binid)  # +1 open end

    # Filter mutations
    header = ['SAMPLE', 'CHROMOSOME', 'POSITION', 'REF', 'ALT', 'BASE', 'CONTEXT', 'SIGNATURE', 'PROB']
    with open(output_f, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        with open(input_f, 'r') as fd:
            next(fd)
            for line in fd:
                # Parse data
                sample, chrom, pos, ref, alt, base, content, sig, prob = line.strip().split('\t')
                for interval in tree[chrom][int(pos)]:
                    ofd.write('{}\n'.format('\t'.join([sample, chrom, pos, ref, alt, base, content, sig, prob])))


if __name__ == '__main__':
    main()

