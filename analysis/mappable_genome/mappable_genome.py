"""Compute mappable genome coordinates"""

import os, sys
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

import click
import pandas as pd
import pybedtools

@click.command()
@click.option('-m', '--mappable_regions_f', default=None, help='Mappable regions')
@click.option('-b', '--mappable_blacklist_f', default=None, help='Blacklisted regions of low mappability')
@click.option('-pv', '--pop_variants_f', default=None, help='Population variants (blacklisted')
@click.option('-o', '--output_file', default=None, help='Output file')
def main(mappable_regions_f, mappable_blacklist_f, pop_variants_f, output_file):
    """
    Create mappable genome annotations: read mappable regions, remove blacklisted regions and population variants

    """

    # Load high mappability regions
    header = ['CHR', 'START', 'END']
    df = pd.read_csv(mappable_regions_f, sep='\t', header=None)
    df.columns = header
    # Transform to BED format (subtract 1 from start, BED is 0 based)
    df['START'] = df.apply(lambda x: x['START'] - 1, axis=1)
    bed_mappable = pybedtools.BedTool.from_dataframe(df)

    # Load blacklisted regions
    blacklist_map_df = pd.read_csv(mappable_blacklist_f, sep='\t', header=None)
    blacklist_map_df.columns = header

    # Load population variants
    blacklist_snp_df = pd.read_csv(pop_variants_f, sep='\t', header=None)
    blacklist_snp_df = blacklist_snp_df[[0,1]]
    blacklist_snp_df[2] = blacklist_snp_df[1]
    blacklist_snp_df.columns = header

    # Concat annotations to exclude
    blacklist_df = pd.concat([blacklist_map_df, blacklist_snp_df])
    blacklist_df.sort_values(by=header, axis=0, ascending=True, inplace=True)
    # Transform to BED format (subtract 1 from start, BED is 0 based)
    blacklist_df['START'] = blacklist_df.apply(lambda x: x['START'] - 1, axis=1)
    bed_blacklist = pybedtools.BedTool.from_dataframe(blacklist_df)
    # Merge
    bed_blacklist_merged = bed_blacklist.merge()

    # Compute mappable genome
    # Subtract blacklisted regions from mappable regions
    mappable_genome_bed = bed_mappable.subtract(bed_blacklist_merged)

    # Return to df and 1-based format
    mappable_genome_df = pd.read_csv(mappable_genome_bed.fn, sep='\t', header=None)
    mappable_genome_df.sort_values(by=[0, 1, 2], axis=0, ascending=True, inplace=True)
    mappable_genome_df.columns = ['CHR', 'START', 'END']
    mappable_genome_df['START'] = mappable_genome_df.apply(lambda x: x['START'] + 1, axis=1)    # 1 based format

    # Save
    mappable_genome_df.to_csv(output_file, sep="\t", compression='gzip', index=None, header=True)


if __name__ == '__main__':
    main()
