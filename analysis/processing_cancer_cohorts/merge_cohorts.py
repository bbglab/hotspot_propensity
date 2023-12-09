"""Merge somatic mutations from individual cohorts into a cancer type"""

import os
import gzip

from bgparsers import readers
import click
import pandas as pd

chromosomes = list(map(str, range(1, 23))) + ['X', 'Y']
genome = 'hg38'
nucleotides = 'ACGT'

@click.command()
@click.option('-i', '--input-directory', default=None, help='Input directory where cohorts are')
@click.option('-a', '--cohorts-annotations', default=None, help='File containing cohort data')
@click.option('-ct', '--cancertype', default=None,
              help='Cancer type name')
@click.option('-o', '--output-file', default=None, help='Output file')
@click.option('-co', '--cohort', default=None, multiple=True, type=click.STRING, help='Cohort to merge')
def main(input_directory, cohorts_annotations, cancertype, output_file, cohort):
    """Merge cohorts into cancer types"""

    cohorts = set(cohort)

    cohorts_annotations_df = pd.read_csv(cohorts_annotations, sep='\t', header=0)
    files_type = input_directory.split('/')[-1]

    header = [
        'CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE',
        'COHORT', 'CANCER_TYPE', 'PLATFORM', 'TYPE', 'AGE', 'TREATED', 'MUTYPE'
    ]
    with gzip.open(output_file, 'wt') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))

        for cohort in cohorts:
            cohort_info_df = cohorts_annotations_df.loc[cohorts_annotations_df['COHORT'] == cohort]
            platform = list(cohort_info_df['PLATFORM'])[0]
            ttype = list(cohort_info_df['TYPE'])[0]
            age = list(cohort_info_df['AGE'])[0]
            treated = list(cohort_info_df['TREATED'])[0]

            if files_type == 'cohorts_raw':
                if 'HARTWIG' in cohort:
                    sufix = '.in.tsv.gz'
                else:
                    sufix = '.in.gz'
            elif files_type == 'cohorts_filtered':
                sufix = '.filtered.in.gz'

            cohort_mutations_file = os.path.join(input_directory, f'{cohort}{sufix}')
            for row in readers.variants(
                    file=cohort_mutations_file,
                    required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE'],
                    extra=['MUTYPE']
            ):
                chrom = row['CHROMOSOME']
                pos = row['POSITION']
                ref = row['REF']
                alt = row['ALT']
                mutype = row['MUTYPE']
                sample = row['SAMPLE']


                ofd.write('{}\n'.format('\t'.join(
                            list(map(str, [chrom, pos, ref, alt, sample,
                                           cohort, cancertype, platform, ttype, age, treated, mutype])))))


if __name__ == '__main__':
    main()
