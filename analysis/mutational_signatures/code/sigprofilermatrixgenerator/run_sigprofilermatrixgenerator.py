"""Run SigProfilerMatrixGenerator (https://github.com/AlexandrovLab/SigProfilerMatrixGenerator)"""

import click

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

@click.command()
@click.option('-i', '--input-directory', default=None, required=True)
def run(input_directory):
    """
     Run SigProfilerMatrixGenerator

     Args:
         input_directory: path to input directory containing a unique mutations file named project.txt

    """

    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        project=input_directory.split('/')[-1],
        genome='GRCh38',
        vcfFiles=input_directory,
        plot=False,
        seqInfo=True
    )


if __name__ == '__main__':
    run()
