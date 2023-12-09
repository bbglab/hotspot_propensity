"""Compute likelihood of a mutation arising from a signature"""

import os
from collections import defaultdict

import click


@click.command()
@click.option('-p', '--context_probs_f', default=None, type=click.Path(exists=True),
              help='Signature probabilities per sample and context')
@click.option('-m', '--muts_dir', default=None, help='Directory containing input mutations')
@click.option('-st', '--sigstype', default=None, help='SBS96 or ID83')
@click.option('-o', '--output_f', default=None, help='Output file')
def main(context_probs_f, muts_dir, sigstype, output_f):
    """"""

    # Load signature probabilities per sample and trinucleotide context
    probs_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    with open(context_probs_f, 'r') as fd:
        for i, line in enumerate(fd):
            # Parse data
            if i > 0:
                data = line.strip().split('\t')
                sample = data[0]
                context = data[1]
                probabilities = data[2:]
                # Load
                for s, p in list(zip(signatures, probabilities)):
                    probs_dict[sample][context][s] = float(p)
            # Parse header
            else:
                header = line.strip().split('\t')
                signatures = header[2:]

    # Read mutations and assign probabilities
    # Write to output file
    if sigstype == 'SBS96':
        header = ['SAMPLE', 'CHROMOSOME', 'POSITION', 'CONTEXT'] + sorted(signatures)
        with open(output_f, 'w') as ofd:
            ofd.write('{}\n'.format('\t'.join(header)))
            for file in os.scandir(muts_dir):
                if file.name.endswith('seqinfo.txt'):
                    with open(file, 'r') as fd:
                        for line in fd:
                            sample, chrom, pos, context, _ = line.strip().split('\t')
                            context_tri = context[3:-1]
                            line = [sample, chrom, pos, context_tri]
                            for signature in sorted(probs_dict[sample][context_tri].keys()):
                                line += [str(probs_dict[sample][context_tri][signature])]
                            ofd.write('{}\n'.format('\t'.join(line)))


if __name__ == '__main__':
    main()

