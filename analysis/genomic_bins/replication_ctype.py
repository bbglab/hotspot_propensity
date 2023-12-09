from collections import defaultdict
import os
import gzip

import click
import pyBigWig
import numpy as np
import pandas as pd


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

cell_lines = {
    'SOLID': ['Helas3','Hepg2','Huvec','Imr90','Mcf7','Nhek','Sknsh'],
    'NON_SOLID': ['Gm06990','Gm12801','Gm12812','Gm12813']
}
CHR = [str(i) for i in range(1, 23)] + ['X', 'Y']

@click.command()
@click.option('--cancer-type', default=None, required=True, type=str)
@click.option('--bins_f', default=None, required=True, type=str)
@click.option('--url', default=None, required=True, type=str)
@click.option('--output_f', default=None, required=True, type=str)
def main(cancer_type, bins_f, url, output_f):

    """Get Replication timing signal from available Roadmap/ENCODE epigenomes"""

    # Output files
    output_log = output_f + '.log'

    # Input bins
    bins_df = pd.read_csv(bins_f, sep='\t', header=0, low_memory=False)
    bins_df.sort_values(by=['CHR', 'START'], ascending=[True, True], inplace=True)

    results_per_bin = defaultdict(list)
    failed_bins = []
    for cell in cell_lines[cancer_type]:
        file = os.path.join(url, f'wgEncodeUwRepliSeq{cell}WaveSignalRep1.bigWig')
        bw = pyBigWig.open(file)

        # Iterate through bins and get data
        for _, x in bins_df.iterrows():
            chrom, start, end, binid = x["CHR"], x["START"], x["END"], x['BINID']
            try:
                bin_data = bw.stats(x["CHR"], x["START"], x["END"], type='mean')
                if bin_data[0]:
                    results_per_bin[binid].append(bin_data[0])
                else:
                    results_per_bin[binid].append(0.0)
                    failed_bins.append((cell, binid))
            except RuntimeError:
                failed_bins.append((cell, binid))

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
