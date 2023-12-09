"""Remove mutations overlapping genomic regions of cancer driver genes"""

from collections import defaultdict
import gzip

from bgparsers import readers
import click
from intervaltree import IntervalTree

@click.command()
@click.option('-m', '--mutations_f', default=None, help='Input file containing somatic mutation of a cancer type')
@click.option('-d', '--drivers_f', default=None, help='Driver genes regions')
@click.option('-o', '--output_f', default=None, help='Output file')
def main(mutations_f, drivers_f, output_f):
    """
    Remove mutations in cancer drivers
    """

    # Load driver regions
    tree = defaultdict(IntervalTree)
    with open(drivers_f, 'r') as fd:
        next(fd)
        for line in fd:
            chrom, start, end = line.strip().split('\t')
            tree[chrom].addi(int(start), int(end) + 1)  # +1 open end
    print('Genomic regions loaded')

    header = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'COHORT', 'CANCER_TYPE', 'PLATFORM', 'TYPE', 'AGE',
              'TREATED', 'MUTYPE']
    extra_cols = ['REF', 'ALT', 'SAMPLE', 'COHORT', 'CANCER_TYPE', 'PLATFORM', 'TYPE', 'AGE', 'TREATED', 'MUTYPE']
    with gzip.open(output_f, 'wt') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))

        for row in readers.variants(
                file=mutations_f,
                required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE'],
                extra=['COHORT', 'CANCER_TYPE', 'PLATFORM', 'TYPE', 'AGE', 'TREATED', 'MUTYPE']
        ):
            chrom = row['CHROMOSOME']
            pos = row['POSITION']

            if tree[chrom][int(pos)]:
                continue
            else:
                ofd.write('{}\n'.format('\t'.join([chrom, str(pos)] + [row[c] for c in extra_cols])))


if __name__ == '__main__':
    main()
