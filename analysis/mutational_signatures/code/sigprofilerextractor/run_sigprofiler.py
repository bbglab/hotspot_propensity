"""
Run Sigprofiler
Full info at: https://osf.io/t6j7u/wiki/3.%20Using%20the%20Tool%20-%20Input/
"""

import click
from SigProfilerExtractor import sigpro as sig

GENOME = 'GRCh38'

@click.command()
@click.option('-i', '--input', default=None, required=True)
@click.option('-o', '--output', default=None, required=True)
@click.option('-c', '--cores', default=None, required=True)
@click.option('-ct', '--context-type', default=None, required=True)
@click.option('--nmf-replicates', default=500, required=True)
@click.option('--max-sigs', default=20, required=True)
def main(input, output, cores, context_type, nmf_replicates, max_sigs):

    sig.sigProfilerExtractor(
        input_type='matrix',
        output=output,
        input_data=input,
        reference_genome=GENOME,
        opportunity_genome=GENOME,
        context_type=context_type,
        exome=False,
        minimum_signatures=1,
        maximum_signatures=max_sigs,
        nmf_replicates=nmf_replicates,
        resample=True,
        seeds='random',
        min_nmf_iterations=10000, max_nmf_iterations=1000000, nmf_test_conv=10000,
        stability=0.8,
        cosmic_version=3.2,
        get_all_signature_matrices=True,
        refit_denovo_signatures=True,
        cpu=int(cores)
    )


if __name__ == '__main__':
    main()
