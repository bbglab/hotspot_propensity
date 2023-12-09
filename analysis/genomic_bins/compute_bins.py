import os
import gzip

import click
from bgreference import refseq
import numpy as np

CHR = [str(i) for i in range(1, 23)] + ['X', 'Y']


def compute_sizes(genome):
    """Given a reference genome, compute the length per chromosome"""
    sizes = {}
    for chr_ in CHR:
        seq = refseq(genome, chr_, start=(1), size=None)
        sizes[chr_] = (1, len(seq))
    return sizes

@click.command()
@click.option('-g', '--genome', default=None, required=True)
@click.option('-bin', '--bin_size', default=None, required=True, type=int)
@click.option('-o', '--output_f', default=None, required=True)
def main(genome, bin_size, output_f):
    """Split the genome in bins (chunks) of equal size"""

    # Compute chromosome size
    chr_sizes = compute_sizes(genome)

    # Split into bin sizes
    header = ['CHR', 'START', 'END']
    with gzip.open(output_f, 'wt') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        for chromosome in CHR:
            start, length = chr_sizes[chromosome]
            start = start - 1  # BED format
            end = length
            positions = np.arange(start, end, bin_size)    # split chromosome into bins
            for i, coordinate in enumerate(positions[1:], 1):
                # Write one bin per line (e.g., "1"\t"0"\t"1000000")
                bin_start = positions[i - 1]
                bin_end = coordinate
                ofd.write('{}\n'.format('\t'.join([chromosome, str(int(bin_start)), str(int(bin_end))])))


if __name__ == '__main__':
    main()
