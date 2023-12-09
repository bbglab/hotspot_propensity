"""Reformat somatic mutation files to fit sigprofilermatrixgenerator input"""

import gzip

from bgparsers import readers
import click


@click.command()
@click.option('-i', '--mutations-file', default=None, required=True)
@click.option('-o', '--output-file', default=None, required=True)
def main(mutations_file, output_file):
    """Reformat file. Include only SNVs"""

    header = [
        'Project',
        'Sample',
        'ID',
        'Genome',
        'mut_type',
        'chrom',
        'pos_start',
        'pos_end',
        'ref',
        'alt',
        'Type'
    ]

    with open(output_file, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))

        with gzip.open(mutations_file, 'rt') as fd:
            project = mutations_file.split('/')[-1].split('.')[0]
            for row in readers.variants(
                    file=mutations_file,
                    required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']
            ):
                chrom = row['CHROMOSOME']
                pos = row['POSITION']
                ref = row['REF']
                alt = row['ALT']
                alt_type = row['ALT_TYPE']
                sample = row['SAMPLE']

                if alt_type == 'snp':  
                    muttype = 'SNP'
                    ofd.write('{}\n'.format('\t'.join(
                        list(map(str, [
                            project, sample, '.', 'GRCh38', muttype, chrom, pos, pos, ref, alt, 'SOMATIC'])))))


if __name__ == '__main__':
    main()

