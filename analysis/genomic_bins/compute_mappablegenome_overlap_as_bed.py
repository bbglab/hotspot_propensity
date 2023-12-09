import os, sys
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

import click
import pandas as pd
import pybedtools

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command()
@click.option('-chr', '--chromosome', default=None, required=True)
@click.option('-b', '--bins_f', default=None, required=True)
@click.option('-g', '--mappable_genome_f', default=None, required=True)
@click.option('-o', '--output_f', default=None, required=True)
def main(chromosome, bins_f, mappable_genome_f, output_f):
    """For a given chromosome and a bin, compute the fraction of the bin that overlaps with the mappable genome"""

    # Read mappable genome coordinates and transform to BED
    mappable_genome_df = pd.read_csv(mappable_genome_f, sep='\t', header=0, low_memory=False)
    mappable_genome_df['START'] = mappable_genome_df.apply(lambda x: x['START'] - 1, axis=1)    # BED format
    mappable_genome_df = mappable_genome_df.loc[mappable_genome_df['CHR'] == chromosome].copy()
    mappable_genome_bed = pybedtools.BedTool.from_dataframe(mappable_genome_df)

    # Read bin coordinates, subset the chromosome under analysis and intersect with the mappable genome
    bins_df = pd.read_csv(bins_f, sep='\t', header=0)
    bins_df['CHR'] = bins_df.apply(lambda x: 'chr' + str(x['CHR']), axis=1)
    bins_df = bins_df.loc[bins_df['CHR'] == chromosome].copy()
    bins_df['BINID'] = bins_df.apply(lambda x: f"{x['CHR']}:{x['START']}-{x['END']}", axis=1)    # bin identifiers
    # Compute overlap for each bin in the chromosome
    with open(output_f, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(['CHR', 'START', 'END', 'BINID', 'BP_OVERLAP'])))
        for binid, data in bins_df.groupby('BINID'):
            bin_bed = pybedtools.BedTool.from_dataframe(data)
            map_bin_bed = bin_bed.intersect(mappable_genome_bed)    # overlap with mappable genome
            # If there is overlap, write the number of bp overlapping the mappable genome
            if len(map_bin_bed) > 0:
                map_bin_df = pd.read_csv(map_bin_bed.fn, sep='\t', header=None)
                map_bin_df.columns = ['CHR', 'START', 'END', 'BINID']
                map_bin_df['LENGTH'] = map_bin_df.apply(lambda x: x['END'] - x['START'], axis=1)    # This is BED format
                total_length = sum(map_bin_df['LENGTH'])
                info = [data['CHR'].iloc[0], str(data['START'].iloc[0]), str(data['END'].iloc[0]), binid, str(total_length)]
            else:
                info = [data['CHR'].iloc[0], str(data['START'].iloc[0]), str(data['END'].iloc[0]), binid, '0']
            ofd.write('{}\n'.format('\t'.join(info)))


if __name__ == '__main__':
    main()
