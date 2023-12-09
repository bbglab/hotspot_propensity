"""Split mutations as inside and outside hotspots"""

import os
import gzip

import click
import pandas as pd

@click.command()
@click.option('-c', '--cancertype', default=None, help='Cancer type name')
@click.option('-m', '--mutations', default=None, type=click.Path(exists=True),
              help='User input file containing somatic mutations in TSV format')
@click.option('-h', '--hotspots', default=None, help='Input hotspots file')
@click.option('-o', '--output-directory', default=None, help='Output file')
def main(cancertype, mutations, hotspots, output_directory):
    """Split mutations into in and out of hotspots"""

    # Outputs
    inhotspot_f = os.path.join(output_directory, f'{cancertype}.mutations_in_hotspots.gz')
    outhotspot_f = os.path.join(output_directory, f'{cancertype}.mutations_out_hotspots.gz')

    # Headers
    header = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'MUTYPE']

    # Write IN mutations
    hotspots_df = pd.read_csv(hotspots, sep='\t', header=0)

    with gzip.open(inhotspot_f, 'wt') as inh_ofd:
        inh_ofd.write('{}\n'.format('\t'.join(header)))
        for mutype in ['snv', 'mnv', 'ins', 'del']:
            mutype_hotspots_df = hotspots_df.loc[hotspots_df['MUT_TYPE'] == mutype].copy()
            for _, row in mutype_hotspots_df.iterrows():
                chrom = str(row['CHROMOSOME'])
                pos = str(row['POSITION'])
                ref = row['REF']
                total_samples_alts = row['MUTATED_SAMPLES_ALTS']

                if mutype != 'del':
                    for sample_alts in total_samples_alts.split(';'):
                        sample, alts = sample_alts.split('::')
                        for alt in alts.split(','):
                            inh_ofd.write('{}\n'.format('\t'.join([chrom, pos, ref, alt, sample, mutype])))
                else:
                    for sample_dels in total_samples_alts.split(';'):
                        sample, dels = sample_dels.split('::')
                        for deletion in dels.split(','):
                            inh_ofd.write('{}\n'.format('\t'.join([chrom, pos, deletion, '-', sample, mutype])))

    # Load mutations IN
    inhotspot_df = pd.read_csv(inhotspot_f, sep='\t', header=0)
    inhotspot_df['MUTS_IN_H'] = inhotspot_df.apply(
        lambda x: '_'.join([str(x['CHROMOSOME']), str(x['POSITION']), x['REF'], x['ALT'], x['SAMPLE'], x['MUTYPE']]), axis=1)
    in_mutations = set(inhotspot_df['MUTS_IN_H'].to_list())

    # Write OUT mutations
    with gzip.open(outhotspot_f, 'wt') as outh_ofd:
        outh_ofd.write('{}\n'.format('\t'.join(header)))

        with gzip.open(mutations, 'rt') as fd:
            next(fd)
            for line in fd:
                chrom, pos, ref, alt, sample = line.strip().split('\t')[:5]
                mutype = line.strip().split('\t')[-1]
                identifier = '_'.join([chrom, pos, ref, alt, sample, mutype])
                if identifier not in in_mutations:
                    outh_ofd.write('{}\n'.format('\t'.join([chrom, pos, ref, alt, sample, mutype])))


if __name__ == '__main__':
    main()
