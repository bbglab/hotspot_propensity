"""Compute the fraction of shared mutations between samples in a file"""

from collections import defaultdict
import gzip

import click
import pandas as pd

@click.command()
@click.option('-i', '--input-file', default=None, type=click.Path(exists=True),
              help='User input file containing somatic mutations in TSV format')
@click.option('-o', '--output-file', default=None, help='Output file')
def main(input_file, output_file):
    """Compute shared mutations for each pair of samples in a file"""

    # Read mutations
    df = pd.read_csv(input_file, header=0, sep='\t', low_memory=False)

    # Add mutation identifier
    df['ID'] = df.apply(lambda x: '{}_{}_{}_{}'.format(x['CHROMOSOME'], x['POSITION'], x['REF'], x['ALT']), axis=1)

    # Compute shared mutations for each pair of samples
    dic_res = defaultdict(lambda: defaultdict(int))
    for i, data in df.groupby(by='ID'):
        sample_list = data['SAMPLE'].tolist()
        for s in sample_list:
            for s2 in sample_list:
                if s2 != s:
                    dic_res[s][s2] += 1

    # Count total mutations per sample
    # Annotate original cohort
    muts_per_sample = {}
    data_per_sample = {}
    for sample, sample_data in df.groupby('SAMPLE'):
        muts_per_sample[sample] = len(sample_data['ID'])
        data_per_sample[sample] = list(sample_data['COHORT'].unique())[0]

    # Write
    header = ['sample1', 'sample2', 'sample1_mutations', 'sample2_mutations', 'sample1_data', 'sample2_data',
            'shared_mutations', 'sample1_fraction_shared', 'sample2_fraction_shared']
    l = []
    for s1, data in dic_res.items():
        s1_mutations = muts_per_sample[s1]
        s1_data = data_per_sample[s1]
        for s2, shared_mut in data.items():
            s2_mutations = muts_per_sample[s2]
            s2_data = data_per_sample[s2]
            l.append(pd.DataFrame([[
                s1, s2, s1_mutations, s2_mutations, s1_data, s2_data,
                shared_mut, shared_mut / s1_mutations, shared_mut / s2_mutations]]))
    if len(l) > 0: 
        df = pd.concat(l)
        df.columns = header
        df['samples'] = df.apply(lambda x: '::'.join(sorted([x['sample1'], x['sample2']])), axis=1)
        df2 = df.drop_duplicates(subset='samples')
        df2.to_csv(output_file, sep='\t', header=True, index=False, compression='gzip')
    else: 
        with gzip.open(output_file, 'wt') as ofd: 
            ofd.write('{}\n'.format('\t'.join(header)))


if __name__ == '__main__':
    main()

