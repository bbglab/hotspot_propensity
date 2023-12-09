"""Remove mutations that do not pass our filters and/or are complex indels"""

from collections import defaultdict
import gzip

import bgreference as bgref
from bgparsers import readers
import click
from intervaltree import IntervalTree

# Global variables
chromosomes = list(map(str, range(1, 23))) + ['X', 'Y']
genome = 'hg38'
nucleotides = 'ACGT'
mappable_regions = 'hg38_100bp.coverage.regions.gz'
blacklisted_regions = 'ENCFF356LFX.bed.gz'
population_variants = 'gnomad.genomes.r3.0.sites.allchr.af_0.01.tsv.gz'


def load_mappability(file, chr_format='chrN'):
    """
    Load mappability regions into intervaltree
    Args:
        file (path): path to file containing mappability data
        chr_format (str): chromosome format used in file. 

    Returns:
        tree (dict): tree with regions, keys are chromosomes, data are regions
    """

    tree = defaultdict(IntervalTree)
    trim = 3 if chr_format == 'chrN' else None
    with gzip.open(file, 'rt') as fd:
        for line in fd:
            chrom, start, end = line.strip().split('\t')
            tree[chrom[trim:]].addi(int(start), int(end) + 1)  # +1 interval
    return tree


def load_variation(file, chr_format='chrN'):
    """
    Load population variants into set
    Args:
        file (path): path to file containing population variants data
        chr_format (str): chromosome format used in file. 

    Returns:
        set_of_interest (set): set with variants annotated as 'N_position'
    """

    set_of_interest = set()
    trim = 3 if chr_format == 'chrN' else None
    with gzip.open(file, 'rt') as fd:
        for line in fd:
            chrom, position = line.strip().split('\t')[:2]
            set_of_interest.add(f'{chrom[trim:]}_{position}')
    return set_of_interest


@click.command()
@click.option('-i', '--input-file', default=None, type=click.Path(exists=True),
              help='User input file containing somatic mutations in TSV format')
@click.option('-o', '--output-file', default=None, help='Output file')
def main(input_file, output_file):
    """Filter a cohort of mutations using HotspotFinder filters"""

    # Load data
    # Mappability data
    mappable_regions_tree = load_mappability(mappable_regions)
    blacklisted_regions_tree = load_mappability(blacklisted_regions)

    # Variation data
    variation_data_set = load_variation(population_variants)

    # Read and filter cohort mutations
    total_mut = 0
    pass_mut = 0
    fail_mut = 0
    header = ['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'MUTYPE']
    with gzip.open(output_file, 'wt') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        for row in readers.variants(
                file=input_file,
                required=['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE']
        ):
            chrom = row['CHROMOSOME']
            pos = row['POSITION']
            chr_pos = f'{chrom}_{pos}'
            ref = row['REF']
            alt = row['ALT']
            alt_type = row['ALT_TYPE']
            sample = row['SAMPLE']
            total_mut += 1

            # Mutations parsed in the same way as HotspotFinder
            # Read mutations in autosomal + sexual chromosomes
            # Mutations that ref != alt are kept
            # Mutations that ref == bgreference are kept
            # Mutations that don't have N in alt/ref nor 5-mer context are kept
            # Complex indels removed

            if chrom in set(chromosomes) and ref != alt:
                fail_filters = 0

                # Read mutations
                # Read substitutions
                if alt_type == 'snp':
                    muttype = 'snv'
                    # Reference filter
                    if ref != bgref.refseq(genome, chrom, pos, 1):
                        fail_filters += 1
                    # Alternate filter
                    if alt not in nucleotides:
                        fail_filters += 1
                    # Context filter
                    pentamer_sequence = bgref.refseq(genome, chrom, pos - 2, 5)
                    if len(pentamer_sequence) != len(
                            [i for i in pentamer_sequence if i in nucleotides]):
                        fail_filters += 1

                # Read MNVs
                elif alt_type == 'mnp':
                    muttype = 'mnv'
                    # Reference filter
                    if ref != bgref.refseq(genome, chrom, pos, len(ref)):
                        fail_filters += 1
                    # Alternate filter
                    if len(alt) != len([i for i in alt if i in nucleotides]):
                        fail_filters += 1
                    # Context filter
                    pentamer_sequence = bgref.refseq(genome, chrom, pos - 2, 5)
                    if len(pentamer_sequence) != len(
                            [i for i in pentamer_sequence if i in nucleotides]):
                        fail_filters += 1

                # Read indels
                elif alt_type == 'indel':

                    # Simple insertions (e.g, G>GA, G>GAAA)
                    if ref == '-':
                        muttype = 'ins'
                        # Alternate filter
                        if len(alt) != len([i for i in alt if i in nucleotides]):
                            fail_filters += 1
                        # Context filter
                        pentamer_sequence = bgref.refseq(genome, chrom, pos - 2, 5)
                        if len(pentamer_sequence) != len(
                                [i for i in pentamer_sequence if i in nucleotides]):
                            fail_filters += 1

                    # Simple deletions (e.g, GT>G, GTG>G)
                    elif alt == '-':
                        muttype = 'del'
                        # Reference filter
                        if ref != bgref.refseq(genome, chrom, pos, len(ref)):
                            fail_filters += 1
                            # Alternate filter
                        if len(ref) != len([i for i in ref if i in nucleotides]):
                            fail_filters += 1
                        # Context filter
                        pentamer_sequence = bgref.refseq(genome, chrom, pos - 2, 5)
                        if len(pentamer_sequence) != len(
                                [i for i in pentamer_sequence if i in nucleotides]):
                            fail_filters += 1
                    # Complex indels
                    else:
                        fail_filters += 1

                # Intersect mappability
                is_mappable = False
                for _ in mappable_regions_tree[chrom][int(pos)]:
                    is_mappable = True
                    break
                for _ in blacklisted_regions_tree[chrom][int(pos)]:
                    fail_filters += 1
                    break

                # Intersect population variants
                if {chr_pos}.intersection(variation_data_set):
                    fail_filters += 1

                # Write mutation if all filters are pass
                if is_mappable and fail_filters == 0:
                    pass_mut += 1
                    ofd.write('{}\n'.format('\t'.join(list(map(str, [chrom, pos, ref, alt, sample, muttype])))))
                else:
                    fail_mut += 1
            else:
                fail_mut += 1

    with open(output_file + '.log', 'w') as ofd:
        for name, info in [('TOTAL', total_mut), ('PASS', pass_mut), ('FAIL', fail_mut)]:
            ofd.write(f'{name}\t{info}\n')


if __name__ == '__main__':
    main()

