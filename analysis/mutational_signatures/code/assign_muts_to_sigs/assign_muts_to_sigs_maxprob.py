"""Assign mutation to the signature with maximum probability"""

import click

import bgreference as bgref


base_dir = {
    'A': 'PURINE',
    'C': 'PYRIMIDINE',
    'G': 'PURINE',
    'T': 'PYRIMIDINE',
}

def rev_comp(seq):
    """Compute reverse complementary of a sequence"""
    comp_nucleotides = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    return ''.join(list(map(lambda x: comp_nucleotides[x], seq[::-1])))

@click.command()
@click.option('-m', '--muts_f', default=None, help='Input mutations file')
@click.option('-o', '--output_f', default=None, help='Output file')
def main(muts_f, output_f):
    """Assign a mutation to the signature with maximum probability"""

    # Load mutations per signature
    with open(output_f, 'w') as ofd:
        header = ['SAMPLE', 'CHROMOSOME', 'POSITION', 'REF', 'ALT', 'BASE', 'CONTEXT', 'SIGNATURE', 'PROB']
        ofd.write('{}\n'.format('\t'.join(header)))

        with open(muts_f, 'r') as fd:
            for i, line in enumerate(fd):

                # Parse data
                if i > 0:
                    data = line.strip().split('\t')
                    sample, chr, pos, context = data[0:4]
                    probabilities = list(map(float, data[4:]))

                    # Check reference
                    reference_n = bgref.refseq('hg38', chr, int(pos), 1)
                    if reference_n == context[2]:
                        alternate_n = context[4]
                        base = base_dir[reference_n]
                    else:
                        alternate_n = rev_comp(context[4])
                        base = base_dir[reference_n]

                    # Get max prob signature
                    sigs_probs = sorted(list(zip(signatures, probabilities)), key=lambda x: x[1], reverse=True)
                    max_sig, max_prob = sigs_probs[0]
                    ofd.write('{}\n'.format('\t'.join(list(map(str, [
                        sample, chr, pos, reference_n, alternate_n, base, context, max_sig, max_prob
                    ])))))

                # Parse header
                else:
                    input_header = line.strip().split('\t')
                    signatures = input_header[4:]


if __name__ == '__main__':
    main()

