import os, sys
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]

import click
import pandas as pd
import pybedtools


cpg_contexts = ['A[C>T]G', 'C[C>T]G', 'G[C>T]G', 'T[C>T]G']


@click.command()
@click.option('--chromosome', default=None, required=True, type=str)
@click.option('--muts_f', default=None, required=True, type=str)
@click.option('--cpg_f', default=None, required=True, type=str)
@click.option('--output_file', default=None, required=True, type=str)
def main(chromosome, muts_f, cpg_f, output_file):

    # Load CpG methylation data for the chromosome
    cpg_df = pd.read_csv(cpg_f, sep='\t', header=0)
    # Subset columns to speed up analysis
    cpg_short_df = cpg_df[['CHR', 'POS']].copy()
    # Transform to BED
    cpg_short_df['END'] = cpg_short_df['POS']
    cpg_short_df['START'] = cpg_short_df.apply(lambda x: x['POS'] - 1, axis=1)
    cpg_short_df = cpg_short_df[['CHR', 'START', 'END']]
    cpg_bed = pybedtools.BedTool.from_dataframe(cpg_short_df)

    # Read mutation data
    muts_df = pd.read_csv(muts_f, sep='\t', header=0)
    # Subset SBS1 mutations in CpG sites
    muts_sbs1_df = muts_df.loc[(muts_df['SIGNATURE'] == 'SBS1') & (muts_df['PROB'] >= 0.5)]
    muts_sbs1_df = muts_sbs1_df.loc[muts_sbs1_df['CONTEXT'].isin(cpg_contexts)]
    # Transform to BED
    muts_sbs1_chr_df = muts_sbs1_df.loc[muts_sbs1_df['CHROMOSOME'] == int(chromosome[3:])].copy()
    muts_sbs1_chr_df['END'] = muts_sbs1_chr_df['POSITION']
    muts_sbs1_chr_df['START'] = muts_sbs1_chr_df.apply(lambda x: x['POSITION'] - 1, axis=1)
    muts_sbs1_chr_df = muts_sbs1_chr_df[['CHROMOSOME', 'START', 'END', 'REF', 'ALT', 'BASE', 'CONTEXT']]
    muts_bed = pybedtools.BedTool.from_dataframe(muts_sbs1_chr_df)

    # Intersect both dataframes and keep the intersection
    # Resulting CpGs are mutated
    mutated_cpgs_bed = cpg_bed.intersect(muts_bed)
    # Reformat
    mutated_cpgs_df = pd.read_csv(mutated_cpgs_bed.fn, sep='\t', header=None)
    mutated_cpgs_df = mutated_cpgs_df.drop_duplicates()    # CpGs mutated more than once, keep only 1 entry
    mutated_cpgs_df.columns = ['CHR', 'START', 'POS']
    mutated_cpgs_df['STATUS'] = ['MUTATED'] * len(mutated_cpgs_df)
    mutated_cpgs_df = mutated_cpgs_df[['CHR', 'POS', 'STATUS']]

    # Intersect both dataframes and keep those only present in the CpG table
    # Resulting CpGs are non mutated
    nomut_cpgs_bed = cpg_bed.intersect(muts_bed, v=True)
    # Reformat
    nomut_cpgs_df = pd.read_csv(nomut_cpgs_bed.fn, sep='\t', header=None)
    nomut_cpgs_df.columns = ['CHR', 'START', 'POS']
    nomut_cpgs_df['STATUS'] = ['NO_MUTATED'] * len(nomut_cpgs_df)
    nomut_cpgs_df = nomut_cpgs_df[['CHR', 'POS', 'STATUS']]

    # Concatenate both results
    results = pd.concat([mutated_cpgs_df, nomut_cpgs_df])
    results.sort_values(by='POS', inplace=True, ascending=True)

    # Append results to original CpG methylation dataframe
    cpg_df.sort_values(by='POS', inplace=True, ascending=True)
    cpg_df = cpg_df.reset_index()
    results = results.reset_index()
    merge = pd.concat([cpg_df, results], axis=1)

    merge.to_csv(output_file, sep='\t', header=True, index=False, compression='gzip')


if __name__ == '__main__':
    main()
