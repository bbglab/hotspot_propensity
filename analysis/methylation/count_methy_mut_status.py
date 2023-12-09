import os, sys
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
import json

import click
import numpy as np
import pandas as pd
import pybedtools

cpg_trinuc = ['ACG', 'CCG', 'GCG', 'TCG']


@click.command()
@click.option('--chromosome', default=None, required=True, type=str)
@click.option('--cpg_f', default=None, required=True, type=str)
@click.option('--bins_trinuc_f', default=None, required=True, type=str)
@click.option('--bins_f', default=None, required=True, type=str)
@click.option('-e', '--epigenomes', default=None, required=True, multiple=True, type=str)
@click.option('--output_file', default=None, required=True, type=str)
def main(chromosome, cpg_f, bins_trinuc_f, bins_f, epigenomes, output_file):

    epigenomes = list(epigenomes)

    # Load CpG methylation data for the chromosome
    cpg_df = pd.read_csv(cpg_f, sep='\t', header=0)
    # Transform to BED
    cpg_df['END'] = cpg_df['POS']
    cpg_df['START'] = cpg_df.apply(lambda x: x['POS'] - 1, axis=1)
    cpg_df = cpg_df[['CHR', 'START', 'END', 'REF_TRI', 'STATUS'] + epigenomes]

    # Remove missing values (annotated as -1.00)
    cpg_filter_df = cpg_df
    for epigenome in epigenomes:
        cpg_filter_df = cpg_filter_df.loc[cpg_filter_df[epigenome] != -1.00].copy()

    # Calculate average methylation fraction per site across epigenomes (always 2 epigenomes)
    cpg_filter_df['ME_MEAN'] = cpg_filter_df.apply(lambda x: np.mean([x[epigenomes[0]], x[epigenomes[1]]]), axis=1)

    # Assign methylation status
    # Methylated cytosines have at least 50% reads methylated
    cpg_filter_df['ME_STATUS'] = cpg_filter_df.apply(lambda x: x['ME_MEAN'] >= 0.5, axis=1)

    # Transform to BED
    cpg_bed = pybedtools.BedTool.from_dataframe(cpg_filter_df)

    # Load dictionary with trinucleotide counts per bin
    with open(bins_trinuc_f, 'r') as fd:
        bins_trinuc = json.load(fd)

    # Load bins
    bins_df = pd.read_csv(bins_f, sep='\t', header=0)
    bins_df['CHR'] = bins_df.apply(lambda x: x['BINID'].split(':')[0], axis=1)
    bins_chr_df = bins_df.loc[bins_df['CHR'] == chromosome]
    bins_chr = bins_chr_df['BINID'].tolist()

    # Create a BED dataframe with bin coordinates
    lines = []
    for binid in bins_chr:
        start, end = binid.split(':')[-1].split('-')
        lines += [pd.DataFrame([[chromosome[3:], start, end, binid]])]
    bins_df = pd.concat(lines)
    bins_df.columns = ['CHR', 'START', 'END', 'BINID']
    bins_bed = pybedtools.BedTool.from_dataframe(bins_df)

    # Get CpGs overlaping bins
    cpgs_bins_bed = cpg_bed.intersect(bins_bed, wao=True)
    cpgs_bins_df = pd.read_csv(cpgs_bins_bed.fn, sep='\t', header=None)
    cpgs_bins_df.columns = list(cpg_filter_df.columns) + list(bins_df.columns) + ['overlap']

    # Save methylation and mutation data per bin
    header = [
        'SIGNATURE',
        'EPIGENOMES',
        'FRAC_METHYLATION',
        'CHR',
        'BINID',
        'TRINUC',
        'TYPE',
        'COUNT'
    ]

    with open(output_file, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        # For each bin
        for binid, binid_cpgs in cpgs_bins_df.groupby('BINID'):
            total_trinuc = bins_trinuc[binid]
            shared_data = ['SBS1', ';'.join(epigenomes), '0.5', chromosome, binid]
            # For each possible CpG trinucleotide
            for trinuc in cpg_trinuc:
                # Total trinucleotide counts
                trinuc_counts = total_trinuc[trinuc]
                trinuc_df = binid_cpgs.loc[binid_cpgs['REF_TRI'] == trinuc]

                methylated_mutated_pos = len(
                    trinuc_df.loc[(trinuc_df['ME_STATUS'] == True) & (trinuc_df['STATUS'] == 'MUTATED')])
                methylated_nomut_pos = len(
                    trinuc_df.loc[(trinuc_df['ME_STATUS'] == True) & (trinuc_df['STATUS'] == 'NO_MUTATED')])
                unmethylated_mutated_pos = len(
                    trinuc_df.loc[(trinuc_df['ME_STATUS'] == False) & (trinuc_df['STATUS'] == 'MUTATED')])
                unmethylated_nomut_pos = len(
                    trinuc_df.loc[(trinuc_df['ME_STATUS'] == False) & (trinuc_df['STATUS'] == 'NO_MUTATED')])

                missing = trinuc_counts - (
                        methylated_mutated_pos + methylated_nomut_pos + unmethylated_mutated_pos + unmethylated_nomut_pos
                )
                # Write
                ofd.write('{}\n'.format('\t'.join(shared_data + [trinuc, 'ME_POS_MUT_POS', str(methylated_mutated_pos)])))
                ofd.write('{}\n'.format('\t'.join(shared_data + [trinuc, 'ME_POS_MUT_NEG', str(methylated_nomut_pos)])))
                ofd.write('{}\n'.format('\t'.join(shared_data + [trinuc, 'ME_NEG_MUT_POS', str(unmethylated_mutated_pos)])))
                ofd.write('{}\n'.format('\t'.join(shared_data + [trinuc, 'ME_NEG_MUT_NEG', str(unmethylated_nomut_pos)])))
                ofd.write('{}\n'.format('\t'.join(shared_data + [trinuc, 'MISSING', str(missing)])))


if __name__ == '__main__':
    main()
