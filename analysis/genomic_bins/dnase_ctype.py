import os
import gzip

import click
import numpy as np
import pandas as pd
import pyBigWig


@click.command()
@click.option('--cancer-type', default=None, required=True, type=str)
@click.option('--bins_f', default=None, required=True, type=str)
@click.option('--metadata_f', default=None, required=True, type=str)
@click.option('--cov_f', default=None, required=True, type=str)
@click.option('--url', default=None, required=True, type=str)
@click.option('--output_f', default=None, required=True, type=str)
def main(cancer_type, bins_f, metadata_f, cov_f, url, output_f):
    """Get DNase signal from available Roadmap/ENCODE epigenomes"""

    # Output files
    output_log = output_f + '.log'

    # Epigenomes
    # Get epigenomes with mark data
    table_df = pd.read_csv(metadata_f, sep='\t', header=0)
    available_epigenomes = list(table_df.loc[table_df['MARK'] == 'DNase']['EID'])

    # Get epigenomes cancer type specific
    table2_df = pd.read_csv(cov_f, sep='\t', header=0)
    ctype_epigenomes = table2_df.loc[table2_df['CANCER_TYPE'] == cancer_type]['DNA_HISTONES_RNA_EPIGENOME'].iloc[0].split(',')

    # Intersect epigenomes for analysis
    epigenomes = list(set(available_epigenomes).intersection(set(ctype_epigenomes)))

    # Input bins
    bins_df = pd.read_csv(bins_f, sep='\t', header=0)
    bins_df.sort_values(by=['CHR', 'START'], ascending=[True, True], inplace=True)

    #### Run
    loaded_epigenomes = []
    for epigenome in epigenomes:
        file = os.path.join(url, f'{epigenome}-DNase.fc.signal.bigwig')
        loaded_epigenomes.append((epigenome, pyBigWig.open(file)))

    header = ['CHR', 'START', 'END', 'ID', 'MEAN_SIGNAL', 'STD']
    failed_bins = []

    with gzip.open(output_f, 'wt') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))

        # Iterate through bins, epigenomes and get data
        for binid, bin_data in bins_df.groupby('BINID'):
            chrom, hg19_start, hg19_end = bin_data["CHR"].iloc[0], bin_data["START"].iloc[0], bin_data["END"].iloc[0]
            hg38_start, hg38_end = binid.split(':')[1].split('-')
            signal = []

            for epi_name, epi_bw in loaded_epigenomes:
                results = epi_bw.stats(chrom, hg19_start, hg19_end, type='mean')
                if results[0]:
                    signal.append(results[0])
                else:
                    signal.append(0.0)
                    failed_bins.append((epi_name, binid))

            # Write
            data_to_write = [
                chrom,
                str(hg38_start),
                str(hg38_end),
                binid,
                str(np.mean(signal)),
                str(np.std(signal))
            ]
            ofd.write('{}\n'.format('\t'.join(data_to_write)))

    with open(output_log, 'w') as ofd:
        for error in failed_bins:
            ofd.write('{}\n'.format('\t'.join(list(error))))


if __name__ == '__main__':
    main()
