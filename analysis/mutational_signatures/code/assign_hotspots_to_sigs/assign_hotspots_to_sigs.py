"""Module to compute the probability of a signature in a hotspot"""

from collections import defaultdict
import gzip

import click


def rev_comp(seq):
    """Compute reverse complementary of a sequence"""
    comp_nucleotides = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    return ''.join(list(map(lambda x: comp_nucleotides[x], seq[::-1])))


sigs_names = {
    'SBS1': 'SBS1',
    'SBS5': 'SBS5',
    'SBS17a': 'SBS17a',
    'SBS17b': 'SBS17b',
    'SBS2': 'APOBEC',
    'SBS13': 'APOBEC',
    'SBS9': 'AID',
    'SBS85': 'AID',
    'SBS84': 'AID',
    'SBS7a': 'UV',
    'SBS7b': 'UV',
    'SBS7c': 'UV',
    'SBS7d': 'UV',
    'SBS4': 'Tobacco',
    'SBS92': 'Tobacco'
}


artifact_sigs = [
    'SBS27',
    'SBS43',
    'SBS45',
    'SBS46',
    'SBS47',
    'SBS48',
    'SBS49',
    'SBS50',
    'SBS51',
    'SBS52',
    'SBS53',
    'SBS54',
    'SBS55',
    'SBS56',
    'SBS57',
    'SBS58',
    'SBS59',
    'SBS60'
]


@click.command()
@click.option('-s', '--context_probs_f', default=None, type=click.Path(exists=True),
              help='Signature probabilities per sample and context')
@click.option('-h', '--hotspots_f', default=None, help='Directory containing input mutations')
@click.option('-o', '--output_f', default=None, help='Output file')
@click.option('-p', '--p_mode', default=None, help='Probabilities as sum or normalized sum')
def main(context_probs_f, hotspots_f, output_f, p_mode):
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
                context_for = context[0] + context[2] + context[-1] + '>' + context[4]
                context_rev = rev_comp(context[0] + context[2] + context[-1]) + '>' + rev_comp(context[4])
                probabilities = data[2:]
                # Load
                for s, p in list(zip(signatures, probabilities)):
                    probs_dict[sample][context_for][s] = float(p)
                    probs_dict[sample][context_rev][s] = float(p)
            # Parse header
            else:
                header = line.strip().split('\t')
                signatures = header[2:]
    noartifact_signatures = list(set(signatures).difference(set(artifact_sigs)))

    # Write to output file
    header = ['HOTSPOT_ID', 'CHROMOSOME', 'POSITION', 'MUT_SAMPLES', 'CONTEXT', 'SIGNATURE', 'NAME', 'PROBABILITY']
    with open(output_f, 'w') as ofd:
        ofd.write('{}\n'.format('\t'.join(header)))

        with gzip.open(hotspots_f, 'rt') as fd:
            next(fd)
            for line in fd:
                # Read only SNVs
                if line.strip().split('\t')[4] == 'snv':
                    chrom, pos, _, hotspot_id = line.strip().split('\t')[:4]
                    mutated_samples = line.strip().split('\t')[21].split(';')
                    context_tri = line.strip().split('\t')[13] + '>' + hotspot_id[-1]

                    # Get probabilities
                    probs_hotspot = defaultdict(float)
                    for sample in mutated_samples:
                        for signature in noartifact_signatures:
                            probs_hotspot[signature] += float(probs_dict[sample][context_tri][signature])
                    # Normalize probabilities if specified
                    if p_mode == 'normsum':
                        total_p = sum(probs_hotspot.values())
                        probs_hotspot = dict([(k, v / total_p) for k, v in probs_hotspot.items()])
                    # Write output
                    for signature in noartifact_signatures:
                        line = list(map(str, [
                            hotspot_id, chrom, pos, len(mutated_samples), context_tri,
                            signature, sigs_names.get(signature, 'Other'), probs_hotspot[signature]
                        ]))
                        ofd.write('{}\n'.format('\t'.join(line)))


if __name__ == '__main__':
    main()
