"""Intersect a given genomic feature and its flanks with mutations"""

from collections import defaultdict

import click
from intervaltree import IntervalTree


@click.command()
@click.option('-m', '--mutations_f', default=None, required=True)
@click.option('-f', '--feature_file', default=None, required=True)
@click.option('-b', '--bp_analysis', default=None, required=True, type=int)
@click.option('-o', '--output_f', default=None, required=True)
def main(mutations_f, feature_file, bp_analysis, output_f):

    #### Load
    tree = defaultdict(IntervalTree)
    with open(feature_file, 'r') as fd:
        next(fd)
        for line in fd:
            chrom, start, end, ident, length = line.strip().split('\t')
            flank_length = (bp_analysis - int(length)) // 2
            right_margin = int(start) - flank_length
            left_margin = int(end) + flank_length
            ident = f'{ident}__{chrom}:{right_margin}-{left_margin}'
            tree[chrom].addi(right_margin, left_margin + 1, ident)  # +1 open end

    #### Mutations
    header = ['ID', 'POS', 'POS_REL_START']

    with open(output_f, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(header + ['SIGNATURE'])))
        with open(mutations_f, 'r') as fd:
            next(fd)
            for line in fd:
                # Parse data
                _, chrom, pos, _, _, _, _, sig, _ = line.strip().split('\t')
                for interval in tree[chrom][int(pos)]:
                    element = interval.data
                    pos_rel = int(pos) - int(element.split('__')[-1].split(':')[-1].split('-')[0])
                    ofd.write('{}\n'.format('\t'.join(list(map(str, [element, pos, pos_rel, sig])))))


if __name__ == '__main__':
    main()
