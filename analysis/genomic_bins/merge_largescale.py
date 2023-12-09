from collections import defaultdict
import gzip

from bgparsers import readers
import click
from intervaltree import IntervalTree
import numpy as np
import pandas as pd


autosomes = list(map(lambda x: f'chr{x}', range(1, 23)))

covariates_f = 'covariates_list.csv'     # this is equivalent to Extended View Dataset EV6 
covariates_df = pd.read_csv(covariates_f, sep='\t', header=0)
cancer_class = {}
for _, row in covariates_df.iterrows():
    if 'Helas3' in row['REPLICATION_EPIGENOME']:
        cancer_class[row['CANCER_TYPE']] = 'SOLID'
    else:
        cancer_class[row['CANCER_TYPE']] = 'NON_SOLID'

@click.command()
@click.option('-c', '--ctype', default=None, required=True)
@click.option('-m', '--mtype', default=None, required=True)
@click.option('-b', '--binsize', default=None, required=True)
@click.option('-o', '--output_path', default=None, required=True)
def main(ctype, mtype, binsize, output_path):

    #### Bins
    bins_dir = '/bins_genome/data'
    bins_file = f'{bins_dir}/hg19_liftoverfromhg38_{binsize}_bin.filtered.bed.gz'
    # Use these bins so that there are no bins with no epigenetic data (data in hg19)

    bins_tree = defaultdict(IntervalTree)
    binids = set()
    with gzip.open(bins_file, 'rt') as fd:
        next(fd)
        for line in fd:
            chrom, _, _, binid = line.strip().split('\t')
            if chrom in autosomes:
                start, end = binid.split(':')[-1].split('-')
                bins_tree[chrom[3:]].addi(int(start), int(end) + 1, binid)  # +1
                binids.add(binid)

    #### DNase
    experiment = 'dnase'
    base_dir = f'/hg38_{binsize}_bins/{experiment}'
    file = f'{base_dir}/{ctype}_hg38_{binsize}_bin.DNase.filtered.bed.gz'
    dnase_df = pd.read_csv(file, sep='\t', header=0)

    #### Histones
    histonemarks = [
        'H3K4me1',
        'H3K4me3',
        'H3K36me3',
        'H3K27me3',
        'H3K9me3',
        'H3K27ac',
        'H3K9ac'
    ]
    experiment = 'histonemarks'
    base_dir = f'/hg38_{binsize}_bins/{experiment}'
    hm_data = {}
    for hm in histonemarks:
        file = f'{base_dir}/{ctype}_hg38_{binsize}_bin.{hm}.filtered.bed.gz'
        df = pd.read_csv(file, sep='\t', header=0)
        hm_data[hm] = df

    #### Replication time
    experiment = 'replication'
    base_dir = f'/hg38_{binsize}_bins/{experiment}'
    metactype = cancer_class[ctype]
    file = f'{base_dir}/{metactype}_hg38_{binsize}_bin.RepliSeq.filtered.bed.gz'
    repli_df = pd.read_csv(file, sep='\t', header=0)

    #### RNA
    experiment = 'rna'
    base_dir = f'/hg38_{binsize}_bins/{experiment}'
    file = f'{base_dir}/{ctype}_hg38_{binsize}_bin.rna.filtered.bed.gz'
    rna_df = pd.read_csv(file, sep='\t', header=0)

    #### Mutations
    base_dir = '/signatures/assign_muts_to_sigs/data'
    sum_probs_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    total_muts = defaultdict(lambda: defaultdict(float))
    if mtype == 'vector':
        signatures = list()
        for sufix in ['in', 'out', 'total']:
            input_f = f'{base_dir}/{ctype}_SBS96_{sufix}.txt'
            with open(input_f, 'r') as fd:
                for i, line in enumerate(fd):
                    # Parse data
                    if i > 0:
                        data = line.strip().split('\t')
                        chrom, pos = data[1:3]
                        probabilities = data[4:]
                        for interval in bins_tree[chrom][int(pos)]:
                            binid = interval.data
                            # Load
                            for s, p in list(zip(signatures, probabilities)):
                                sum_probs_dict[sufix][binid][s] += float(p)
                                total_muts[sufix][binid] += 1
                    # Parse header
                    else:
                        header = line.strip().split('\t')
                        signatures = sorted(header[4:])
    elif mtype == 'maxprob':
        signatures = set()
        for sufix in ['in', 'out', 'total']:
            input_f = f'{base_dir}/{ctype}_SBS96_{sufix}_maxprob.tsv'
            with open(input_f, 'r') as fd:
                next(fd)
                for line in fd:
                    # Parse data
                    _, chrom, pos, _, s, _ = line.strip().split('\t')
                    signatures.add(s)
                    for interval in bins_tree[chrom][int(pos)]:
                        binid = interval.data
                        sum_probs_dict[sufix][binid][s] += 1
                        total_muts[sufix][binid] += 1

    #### Hotspots
    base_dir = '/hotspots/'
    file = f'{base_dir}/{ctype}.results.tsv.gz'

    bins_hotspots = defaultdict(int)
    for row in readers.variants(
            file=file,
            required=['CHROMOSOME', 'POSITION'],
            extra=['MUT_TYPE']

    ):
        chrom = row['CHROMOSOME']
        pos = row['POSITION']
        mutype = row['MUT_TYPE']

        if mutype == 'snv':
            for interval in bins_tree[chrom][int(pos)]:
                bins_hotspots[interval.data] += 1

    print('Data loaded')
    ## Merge
    header = [
                 'DATA_TYPE',
                 'BIN_SIZE',
                 'CHR',
                 'START',
                 'END',
                 'BINID',
                 'CTYPE',
                 'DHS'] + sorted(histonemarks) + [
                 'REPLI_T',
                 'EXPR',
                 'H_SNV_BIN',
                 'M_IN_BIN',
                 'M_OUT_BIN',
                 'M_TOTAL_BIN',
                 'SIG',
                 'M_IN_SIG',
                 'M_OUT_SIG',
                 'M_TOTAL_SIG'
             ]

    lines = []
    for binid in binids:
        binid_data = []
        chrom = binid.split(':')[0]
        start, end = binid.split(':')[-1].split('-')
        binid_data += [mtype, binsize, chrom, start, end, binid, ctype]

        # DNase
        try:
            dnase = round(dnase_df.loc[dnase_df['ID'] == binid]['MEAN_SIGNAL'].iloc[0], 6)
        except IndexError:
            dnase = np.nan
        binid_data += [dnase]

        # Histone marks
        for hm in sorted(histonemarks):
            data = hm_data[hm]
            try:
                binid_data += [round(data.loc[data['ID'] == binid]['MEAN_SIGNAL'].iloc[0], 6)]
            except IndexError:
                binid_data += [np.nan]

        # Replication
        try:
            binid_data += [round(repli_df.loc[repli_df['ID'] == binid]['MEAN_SIGNAL'].iloc[0], 6)]
        except IndexError:
            binid_data += [np.nan]

        # RNA
        try:
            binid_data += [round(rna_df.loc[rna_df['ID'] == binid]['MEAN_SIGNAL'].iloc[0], 6)]
        except IndexError:
            binid_data += [np.nan]

        # Hotspots
        binid_data += [bins_hotspots.get(binid, 0)]

        # Mutations in the bin
        for sufix in ['in', 'out', 'total']:
            total_muts_sufix = total_muts[sufix]
            binid_data += [total_muts_sufix.get(binid, 0)]

        # Mutations in the bin, split by signature
        for sig in signatures:
            sig_data_to_write = [sig]
            for sufix in ['in', 'out', 'total']:
                muts_sufix = sum_probs_dict[sufix]
                muts_sufix_bin = muts_sufix.get(binid, {})
                if muts_sufix_bin:
                    sig_data_to_write += [muts_sufix_bin.get(sig, 0)]
                else:
                    sig_data_to_write += [0]

            lines.append(pd.DataFrame([binid_data + sig_data_to_write]))

    results_df = pd.concat(lines)
    results_df.columns = header

    # Save
    output_f = f'{output_path}/{ctype}_{binsize}_largescale_matrix_{mtype}.tsv'
    results_df.to_csv(output_f, sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()
