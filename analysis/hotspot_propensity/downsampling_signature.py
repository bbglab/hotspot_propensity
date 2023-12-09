"""Subsample a cancer type"""

import os

import click
import numpy as np
np.random.seed(21736512)
import pandas as pd

@click.command()
@click.option('-c', '--cancertype', default=None, help='Cancer type name')
@click.option('-sig', '--signature', default=None, help='Signature')
@click.option('-m', '--mutations', default=None, type=click.Path(exists=True),
              help='User input file containing somatic mutations in TSV format')
@click.option('-o', '--output-directory', default=None, help='Output directory')
@click.option('-i', '--iterations', default=None, help='Number of subsamples')
@click.option('-s', '--sample_size', default=None, help='Number of samples per iteration')
@click.option('-mc', '--mcutoff', default=None, help='Number of mutations selected per sample per iteration')
def main(cancertype, signature, mutations, output_directory, iterations, sample_size, mcutoff):
    """Select n random samples from a cancer type and subset N random mutations"""

    # Create output directory
    output_dir = os.path.join(output_directory, cancertype)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Get subsamples and generate subsamples as individual files containing mutations
    # Get samples with at least 'mcutoff' mutations from signature under analysis
    # Save samples per iteration for the record
    lines = []
    ctype_mutations_df = pd.read_csv(mutations, sep='\t', header=0, low_memory=False)
    ctype_mutations_df = ctype_mutations_df.loc[(ctype_mutations_df['SIGNATURE'] == signature) & (ctype_mutations_df['PROB'] > 0.5)].copy()
    muts_per_sample = dict(ctype_mutations_df['SAMPLE'].value_counts())
    muts_per_sample_filter = dict([(k, v) for k, v in muts_per_sample.items() if v > int(mcutoff)])
    filtered_samples = list(muts_per_sample_filter.keys())
    # Reformat
    ctype_mutations_df = ctype_mutations_df[['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']]
    header = ctype_mutations_df.columns

    for i in range(0, int(iterations)):
        iter_id = f'{cancertype}_iter_{i}'

        # Select a number of samples randomly without replacement (a sample cannot be selected twice in the subsample)
        subsample = np.random.choice(filtered_samples, size=int(sample_size), replace=False)

        # Randomize mutations
        iter_mutations_df = []
        for sample in subsample:
            sample_df = ctype_mutations_df.loc[ctype_mutations_df['SAMPLE'] == sample].copy()
            # Select a number of mutations randomly without replacement (a mutation cannot be selected twice)
            sample_muts_df = sample_df.sample(n=int(mcutoff), replace=False).copy()
            iter_mutations_df.append(sample_muts_df)
        iter_mutations_df = pd.concat(iter_mutations_df)
        iter_mutations_df.to_csv(os.path.join(output_dir, f'{iter_id}_{signature}.in.gz'), sep='\t', header=header, index=0)
        lines += [pd.DataFrame([[f'{iter_id}_{signature}', ','.join(subsample)]])]

    subsampling_df = pd.concat(lines)
    subsampling_df.to_csv(os.path.join(output_dir, f'subsampling_{signature}.txt'), sep='\t', header=0, index=0)
    print('Subsampling generated!')


if __name__ == '__main__':
    main()
