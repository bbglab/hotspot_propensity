import gzip

from bgreference import hg38
import click
import pandas as pd
from pyliftover import LiftOver
lo = LiftOver('hg19', 'hg38')


def rev_comp(seq):
    """Compute reverse complementary of a sequence"""
    comp_nucleotides = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }

    return ''.join(list(map(lambda x: comp_nucleotides[x], seq[::-1])))


autosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']

@click.command()
@click.option('--chrom', default=None, required=True, type=str)
@click.option('--input_file', default=None, required=True, type=str)
@click.option('--output_file', default=None, required=True, type=str)
@click.option('--mnemonics', default=None, required=True, type=str)
def main(chrom, input_file, output_file, mnemonics):
    """
    Liftover CpG sites from hg19 to hg38

    Output file contains successfuly liftovered CpG sites (CpG site in hg38)
    Statistics are saved as [output_file] + .log

    """

    log_output = output_file + '.log'

    columns_df = pd.read_csv(mnemonics, sep='\t', header=None)
    header = ['CHR', 'POS', 'POS_hg19', 'REF', 'REF_TRI'] + columns_df[0].tolist()

    total = 0
    nocpg = 0
    nolift = 0
    success = 0
    with gzip.open(output_file, 'wt') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))
        with open(input_file, 'r') as fd:
            for line in fd:
                total += 1
                pos_hg19 = int(line.strip().split()[0])
                methylation = line.strip().split()[1:]

                # Liftover to hg38
                liftover = lo.convert_coordinate(f'chr{chrom}', pos_hg19)

                if liftover:
                    pos_hg38 = liftover[0][1]

                    # Get sequence
                    try:
                        seq_hg38 = hg38(chrom, pos_hg38 - 1, size=3)
                    except ValueError:
                        nolift += 1
                        continue

                    if 'N' not in seq_hg38:

                        # If CpG in forward
                        if seq_hg38[1] == 'C':
                            if seq_hg38[2] == 'G':
                                data = [chrom, str(pos_hg38), str(pos_hg19), seq_hg38[1], seq_hg38] + methylation
                                ofd.write('{}\n'.format('\t'.join(data)))
                                success += 1
                            else:
                                nocpg += 1

                        # If CpG in reverse
                        elif seq_hg38[1] == 'G':
                            seq_hg38 = rev_comp(seq_hg38)
                            if seq_hg38[2] == 'G':
                                data = [chrom, str(pos_hg38), str(pos_hg19), seq_hg38[1], seq_hg38] + methylation
                                ofd.write('{}\n'.format('\t'.join(data)))
                                success += 1
                            else:
                                nocpg += 1
                        else:
                            nocpg += 1
                    else:
                        nolift += 1
                else:
                    nolift += 1
    
    with open(log_output, 'w') as ofd:
        ofd.write(f'TOTAL\t{total}\n')
        ofd.write(f'PASS\t{success}\n')
        ofd.write(f'NO_LIFTOVER\t{nolift}\n')
        ofd.write(f'NO_CPG\t{nocpg}\n')


if __name__ == '__main__':
    main()
