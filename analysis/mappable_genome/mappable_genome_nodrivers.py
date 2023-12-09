"""Substract driver coordinates from the mappable genome"""

import os, sys
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

import click
import pandas as pd
import pybedtools

@click.command()
@click.option('-m', '--mappable_genome_f', default=None, help='Mappable genome regions')
@click.option('-d', '--drivers_f', default=None, help='Driver gene regions')
@click.option('-o', '--output_file', default=None, help='Output file')
def main(mappable_genome_f, drivers_f, output_file):
    """
    Remove driver coordinate annotations from the mappable genome
    """

    # Load mappable genome
    df = pd.read_csv(mappable_genome_f, sep='\t', header=0)
    # Transform to BED format (subtract 1 from start, BED is 0 based)
    df['START'] = df.apply(lambda x: x['START'] - 1, axis=1)
    bed_mappable = pybedtools.BedTool.from_dataframe(df)

    # Load driver regions
    drivers_df = pd.read_csv(drivers_f, sep='\t', header=0)
    drivers_df['CHR'] = [f'chr{c}' for c in drivers_df['CHR'].tolist()]
    # Transform to BED format (subtract 1 from start, BED is 0 based)
    drivers_df['START'] = drivers_df.apply(lambda x: x['START'] - 1, axis=1)
    drivers_bed = pybedtools.BedTool.from_dataframe(drivers_df)

    # Subtract driver regions from mappable genome
    mappable_genome_bed = bed_mappable.subtract(drivers_bed)

    # Return to df and 1-based format
    mappable_genome_df = pd.read_csv(mappable_genome_bed.fn, sep='\t', header=None)
    mappable_genome_df.sort_values(by=[0, 1, 2], axis=0, ascending=True, inplace=True)
    mappable_genome_df.columns = ['CHR', 'START', 'END']
    mappable_genome_df['START'] = mappable_genome_df.apply(lambda x: x['START'] + 1, axis=1)    # 1 based format

    # Save
    mappable_genome_df.to_csv(output_file, sep="\t", compression='gzip', index=None, header=True)


if __name__ == '__main__':
    main()
