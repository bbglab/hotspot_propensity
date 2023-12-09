"""Compute number of mutations per sample"""

from collections import defaultdict

from bgparsers import readers
import click


@click.command()
@click.option('-i', '--input-file', default=None, required=True)
@click.option('-o', '--output-file', default=None, required=True)
def main(input_file, output_file):
    """Compute the number of mutations per sample"""

    # Compute
    print('Reading mutations...')
    muts_per_sample = defaultdict(lambda: defaultdict(int))
    cohorts_d = defaultdict(set)
    type_d = defaultdict(set)
    age_d = defaultdict(set)
    for row in readers.variants(
            file=input_file,
            required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE'],
            extra=['COHORT', 'TYPE', 'AGE', 'MUTYPE']
    ):
        sample = row['SAMPLE']
        mutype = row['MUTYPE']
        muts_per_sample[sample][mutype] += 1
        cohorts_d[sample].add(row['COHORT'])
        type_d[sample].add(row['TYPE'])
        age_d[sample].add(row['AGE'])

    # Write output
    print('Writing output...')
    header = ['SAMPLE', 'COHORT', 'TYPE', 'AGE', 'SNV', 'MNV', 'INS', 'DEL', 'TOTAL']
    with open(output_file, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        for sample, sample_data in muts_per_sample.items():
            counts = []
            for mutype in ['snv', 'mnv', 'ins', 'del']:
                counts += [sample_data[mutype]]
            ofd.write('{}\n'.format('\t'.join(
                list(map(str, [sample, list(cohorts_d[sample])[0], list(type_d[sample])[0], list(age_d[sample])[0]] + counts + [sum(counts)])))))
    print('Finished')


if __name__ == '__main__':
    main()
