import os, sys
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

import click
import pandas as pd
import pybedtools


@click.command()
@click.option('-b', '--bins_f', default=None, required=True)
@click.option('-g', '--mappable_genome_f', default=None, required=True)
@click.option('-o', '--output_f', default=None, required=True)
def main(bins_f, mappable_genome_f, output_f):
    """
    Intersect the mappable genome with filtered bins to find mappable positions per bin
    """

    # Read mappable genome coordinates and transform to BED
    mappable_genome_df = pd.read_csv(mappable_genome_f, sep='\t', header=0, low_memory=False)
    mappable_genome_df['START'] = mappable_genome_df.apply(lambda x: x['START'] - 1, axis=1)    # BED format
    mappable_genome_bed = pybedtools.BedTool.from_dataframe(mappable_genome_df)

    # Read bin coordinates
    bins_df = pd.read_csv(bins_f, sep='\t', header=0)
    bins_bed = pybedtools.BedTool.from_dataframe(bins_df)

    # Intersect
    intersect_bed = mappable_genome_bed.intersect(bins_bed, wb=True)    # wb writes the original entry in bins_bed
    intersect_df = pd.read_csv(intersect_bed.fn, sep='\t', header=None)

    # Reformat
    intersect_df = intersect_df[[0, 1, 2, 6]]
    intersect_df.columns = ['CHR', 'START', 'END', 'BINID']

    # Save
    intersect_df.to_csv(output_f, sep='\t', index=False)


if __name__ == '__main__':
    main()
