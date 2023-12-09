from collections import defaultdict
import os
import gzip

import click
import numpy as np
import pandas as pd
import pyBigWig


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

available_epigenomes = [
    'E003',
    'E024',
    'E007',
    'E013',
    'E012',
    'E011',
    'E004',
    'E005',
    'E006',
    'E062',
    'E037',
    'E038',
    'E047',
    'E050',
    'E055',
    'E056',
    'E059',
    'E061',
    'E057',
    'E058',
    'E028',
    'E027',
    'E054',
    'E053',
    'E112',
    'E071',
    'E070',
    'E082',
    'E100',
    'E104',
    'E095',
    'E105',
    'E065',
    'E085',
    'E084',
    'E109',
    'E106',
    'E079',
    'E094',
    'E097',
    'E087',
    'E066',
    'E098',
    'E096',
    'E113',
    'E114',
    'E116',
    'E117',
    'E118',
    'E119',
    'E120',
    'E122',
    'E123',
    'E127',
    'E128'
]

no_strand_epigenomes = ['E028', 'E037', 'E038', 'E047', 'E050', 'E062', 'E084', 'E085']

@click.command()
@click.option('--cancer-type', default=None, required=True, type=str)
@click.option('--bins_f', default=None, required=True, type=str)
@click.option('--metadata_f', default=None, required=True, type=str)
@click.option('--cov_f', default=None, required=True, type=str)
@click.option('--url', default=None, required=True, type=str)
@click.option('--output_f', default=None, required=True, type=str)
def main(cancer_type, bins_f, metadata_f, cov_f, url, output_f):
    """Get RNA coverage signal from available Roadmap/ENCODE epigenomes"""

    # Output files
    output_log = output_f + '.log'

    metadata_df = pd.read_csv(metadata_f, sep='\t', header=0)
    epigenome_names = dict(
        zip(metadata_df['Epigenome ID (EID)'], metadata_df['Epigenome name (from EDACC Release 9 directory)']))

    # Get epigenomes cancer type specific
    table2_df = pd.read_csv(cov_f, sep='\t', header=0)
    ctype_epigenomes = table2_df.loc[table2_df['CANCER_TYPE'] == cancer_type]['DNA_HISTONES_RNA_EPIGENOME'].iloc[0].split(',')

    # Intersect epigenomes for analysis
    epigenomes = list(set(available_epigenomes).intersection(set(ctype_epigenomes)))

    # Input bins
    bins_df = pd.read_csv(bins_f, sep='\t', header=0)
    bins_df.sort_values(by=['CHR', 'START'], ascending=[True, True], inplace=True)

    results_per_bin = defaultdict(list)
    failed_bins = []
    for epigenome in epigenomes:
        fullname = epigenome_names[epigenome]

        # Define file and analysis type (stranded or not)
        if epigenome not in set(no_strand_epigenomes):
            rna_pos_file = os.path.join(url, f'{epigenome}.{fullname}.norm.pos.bw')
            rna_neg_file = os.path.join(url, f'{epigenome}.{fullname}.norm.neg.bw')
            rna_pos_bw = pyBigWig.open(rna_pos_file)
            rna_neg_bw = pyBigWig.open(rna_neg_file)

            # Iterate through bins and get data
            for _, x in bins_df.iterrows():
                chrom, start, end, binid = x["CHR"], x["START"], x["END"], x['BINID']

                try:
                    bin_data_pos = rna_pos_bw.stats(chrom, start, end, type='mean')
                    bin_data_neg = rna_neg_bw.stats(chrom, start, end, type='mean')

                    if bin_data_pos[0]:
                        bin_data_pos = bin_data_pos[0]
                    else:
                        bin_data_pos = 0.0
                        failed_bins.append((epigenome, binid, 'pos'))
                    if bin_data_neg[0]:
                        bin_data_neg = -bin_data_neg[0]    # values are given in negative, check README at FTP
                    else:
                        bin_data_neg = 0.0
                        failed_bins.append((epigenome, binid, 'neg'))

                    # Add up pos and neg strands
                    bin_data = bin_data_pos + bin_data_neg
                    results_per_bin[binid].append(bin_data)

                except RuntimeError:
                    failed_bins.append((epigenome, binid, 'pos_neg'))
        else:
            file = os.path.join(url, f'{epigenome}.{fullname}.norm.bw')
            bw = pyBigWig.open(file)

            # Iterate through bins and get data
            for _, x in bins_df.iterrows():
                chrom, start, end, binid = x["CHR"], x["START"], x["END"], x['BINID']
                try:
                    bin_data = bw.stats(chrom, start, end, type='mean')
                    if bin_data[0]:
                        results_per_bin[binid].append(bin_data[0])
                    else:
                        results_per_bin[binid].append(0.0)
                        failed_bins.append((epigenome, binid, 'unknown'))
                except RuntimeError:
                    failed_bins.append((epigenome, binid, 'unknown'))

    header = ['CHR', 'START', 'END', 'ID', 'MEAN_SIGNAL', 'STD']
    with gzip.open(output_f, 'wt') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        for binid, values in results_per_bin.items():
            chromosome = binid.split(':')[0]
            start, end = binid.split(':')[1].split('-')
            print(binid, values)
            data = [
                chromosome,
                start,
                end,
                binid,
                str(np.mean(values)),
                str(np.std(values))
            ]
            ofd.write('{}\n'.format('\t'.join(data)))

    with open(output_log, 'w') as ofd:
        for error in failed_bins:
            ofd.write('{}\n'.format('\t'.join(list(error))))


if __name__ == '__main__':
    main()
