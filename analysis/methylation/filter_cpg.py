import gzip

import click
from collections import defaultdict
from intervaltree import IntervalTree
import pandas as pd

@click.command()
@click.option('--input_file', default=None, required=True, type=str)
@click.option('--output_file', default=None, required=True, type=str)
@click.option('--mappable_f', default=None, required=True, type=str)
def main(input_file, mappable_f, output_file):
    """
    Intersect CpG sites with mappable (driver free) positions
    """

    log_output = output_file + '.log'

    tree = defaultdict(IntervalTree)
    with gzip.open(mappable_f, 'rt') as fd:
        next(fd)  # skip header
        for line in fd:
            chrom, start, end, _ = line.strip().split('\t')
            tree[chrom].addi(int(start) + 1, int(end) + 1)
            
    # Drop duplicates (two hg19 positions mapping to the same hg38 position)
    # Keep data from the first CpG in hg19
    cpg_df = pd.read_csv(input_file, sep='\t', header=0)
    nodups_cpg_df = cpg_df.drop_duplicates(subset=['POS'], keep='first')
    nodups_f = f'{input_file.split(".")[0]}.nodups.txt.gz'
    nodups_cpg_df.to_csv(nodups_f, sep='\t', index=False, header=True, compression='gzip')

    # Intersect CpG with mappable genome
    columns_df = pd.read_csv(input_file, sep='\t', header=0, nrows=1)
    header = columns_df.columns.tolist()
    total = 0
    mappable = 0
    with gzip.open(output_file, 'wt') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        with gzip.open(nodups_f, 'rt') as fd:
            next(fd)    # skip header
            for line in fd:
                total += 1
                chromosome = 'chr' + line.strip().split()[0]
                position = int(line.strip().split()[1])
                if tree[chromosome][position]:
                    mappable += 1
                    ofd.write('{}\n'.format('\t'.join(line.strip().split())))

    with open(log_output, 'w') as ofd:
        ofd.write(f'TOTAL\t{total}\n')
        ofd.write(f'MAPPABLE\t{mappable}\n')

if __name__ == '__main__':
    main()
