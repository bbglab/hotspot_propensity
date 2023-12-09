import os

import click
import pandas as pd


@click.command()
@click.option('-s', '--subsamples', default=None, help='Input directory')
@click.option('-i', '--input_dir', default=None, help='Input directory')
@click.option('-o', '--output_f', default=None, help='Output file')
def main(subsamples, input_dir, output_f):
    """Compute number of hotspots per subsample"""

    subsampling_df = pd.read_csv(subsamples, sep='\t', header=None)

    # Read iterations
    results_table = []
    for i, row in subsampling_df.iterrows():
        iter_name = row[0]

        # Compute hotspots per sample
        iter_samples = row[1].split(',')
        iter_df = pd.read_csv(os.path.join(input_dir, f'{iter_name}.results.tsv.gz'), sep='\t', header=0)
        hotspot_positions_unique = len(iter_df['CHR_POS'].unique())

        results_table += [pd.DataFrame([[
            iter_name,
            len(iter_samples),
            len(iter_df),
            hotspot_positions_unique

        ]])]

    # Save
    results_table = pd.concat(results_table)
    results_table.columns = [
        'ITER',
        'N_SAMPLES',
        'TOTAL_HOTSPOTS',
        'TOTAL_HOTSPOTS_UNIQPOS'
    ]
    results_table.to_csv(output_f, sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()
