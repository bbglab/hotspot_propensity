"""Test the significance of the fold change enrichment/depletion in a chromatin feature"""

from collections import Counter
from collections import defaultdict
import json
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from bgreference import hg38
import click
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV



def rev_comp(seq):
    """Compute reverse complementary of a sequence"""
    comp_nucleotides = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    return ''.join(list(map(lambda x: comp_nucleotides[x], seq[::-1])))


class AnalyticalPvalue:
    """Class to calculate the analytical pvalue of elements"""

    def __init__(self, jobs=1):
        """Initialize the AnalyticalPvalue class

        Args:
            jobs (int): cores to use in the computation

        Returns:
            None
        """

        self.jobs = jobs
        self.gkde = None
        self.bandwidth = None

    def _get_best_estimator(self, expected):
        """Calculate the best bandwidth estimator.

        Args:
            expected: list of floats, list of expected values

        Returns:
            float, best bandwidth estimator
        """

        # parameters of the fitting
        params = {
            'bandwidth': np.logspace(-1, 0, 10),
            'kernel': ['gaussian']
        }
        grid = GridSearchCV(KernelDensity(), params, n_jobs=self.jobs)
        data_newaxis = expected[:, np.newaxis]
        grid.fit(data_newaxis)
        return grid.best_estimator_.bandwidth

    def _get_min_analytical_pvalue(self, up, down, calls=0, pvals=None):
        """Get the best pvalue with a resolution > 0.

        Args:
            up (numeric): higher score
            down (numeric): lower score
            calls (int): number of times the function has been called
            pvals (list): list of pvalues

        Returns:
            float, pvalue > 0
        """

        if pvals is None:
            pvals = []
        calls += 1
        mid = down + ((up - down) / 2)
        analytical_pvalue = self.gkde.integrate_box_1d(mid, float('inf'))
        if analytical_pvalue > 0:
            pvals.append(analytical_pvalue)
        if calls == 100:
            return pvals[-1]
        if analytical_pvalue > 0:
            return self._get_min_analytical_pvalue(up=up, down=mid, calls=calls, pvals=pvals)
        else:
            return self._get_min_analytical_pvalue(up=mid, down=down, calls=calls, pvals=pvals)

    def calculate_bandwidth(self, expected, pseudocount=1):
        """Calculate the best estimator
        Args:
            expected (iterable): array of simulated mutations
            pseudocount (int): value to add when expected are zeroes

        Returns:
            None

        """

        # Add a pseudocount to get rid of LinAlgError("singular matrix")
        if sum(expected) == 0:
            expected = np.array([0] * 999 + [pseudocount])
        else:
            expected = np.array(expected)
        try:
            self.gkde = gaussian_kde(expected)
        except np.linalg.linalg.LinAlgError as e:
            expected[-1] += 10**(-14)
            self.gkde = gaussian_kde(expected)
        self.bandwidth = self._get_best_estimator(expected)
        self.gkde.set_bandwidth(bw_method=self.bandwidth)

    def _calculate_analytical_pvalue(self, observed):
        """Calculate the analytical pvalue of an element

        Args:
           observed (float): observed value

        Returns:
            float, analytical pvalue

        """

        analytical_pvalue = self.gkde.integrate_box_1d(observed, float('inf'))
        if analytical_pvalue == 0:
            analytical_pvalue = self._get_min_analytical_pvalue(up=observed, down=0)
        return analytical_pvalue

    def calculate(self, observed):
        """Calculate the analytical pvalue

        Args:
           observed (float): observed value

        Returns:
            analytical_pvalue (float): analytical pvalue
        """

        analytical_pvalue = self._calculate_analytical_pvalue(observed)
        return analytical_pvalue


def compute_fold_change(data, feature_size, flank_size, total_features):
    """Compute fold change in proportion of mutations in feature vs flank"""

    # Get mutations in the feature and flanks
    data_flank_5 = sum(data[0:flank_size])
    data_feature = sum(data[flank_size:flank_size + feature_size])
    data_flank_3 = sum(data[flank_size + feature_size:])
    data_flank = data_flank_5 + data_flank_3

    # Get number of positions analysed (feature sites)
    population_flank = flank_size * 2 * total_features
    population_feature = feature_size * total_features

    # Calculate fold change OBSERVED
    prop_feature = data_feature / population_feature
    prop_flank = data_flank / population_flank

    return prop_feature / prop_flank


@click.command()
@click.option('-f', '--feature_file', default=None, required=True)
@click.option('-fs', '--feature_size', default=None, required=True, type=int)
@click.option('-ws', '--window_size', default=None, required=True, type=int)
@click.option('-sig', '--signature', default=None, required=True)
@click.option('-sigf', '--signatures_f', default=None, required=True)
@click.option('-nsim', default=None, required=True, type=int)
@click.option('-o', '--output_f', default=None, required=True)
@click.option('--seed', type=click.INT, default=None, help='Seed to use in the simulations')
def main(feature_file, feature_size, window_size, signature, signatures_f, nsim, output_f, seed):

    np.random.seed(seed)
    flank_size = int((window_size - feature_size) / 2)

    # Save reference trinucleotide probabilities (forward and reverse complementary) in a dictionary
    # The probability of a trinucleotide is the sum of the probabilities of the three trinucleotide>alternate
    sbs96_df = pd.read_csv(signatures_f, sep='\t', header=0)
    sbs96_df['Trinucleotide'] = sbs96_df.apply(
        lambda x: x['Type'][0] + x['Type'][2] + x['Type'][-1], axis=1)
    trinuc_probs_dict = defaultdict(float)
    for trinuc, alternates in sbs96_df.groupby('Trinucleotide'):
        trinuc_rev = rev_comp(trinuc)
        sig_sum_probs = alternates[signature].sum()
        trinuc_probs_dict[trinuc] = sig_sum_probs
        trinuc_probs_dict[trinuc_rev] = sig_sum_probs

    # Get observed data and compute observed fold change feature vs flanks
    df = pd.read_csv(feature_file, sep='\t', header=0)
    df = df.loc[df['SIGNATURE'] == signature].copy()
    observed_in = []
    for i in range(0, 2000):
        df_pos = df.loc[df['POS_REL_START'] == i].copy()
        observed_in += [len(df_pos)]
    # Fold change and Chi-square
    total_features = len(df['ID'].unique())
    observed_fc = compute_fold_change(observed_in, feature_size, flank_size, total_features)

    # Simulate mutations
    total_simulations = defaultdict(lambda: defaultdict(int))
    for ge_id, ge_data in df.groupby('ID'):
        ge_mutations = len(ge_data)

        # Genomic element + flanks trinucleotide sequence
        ge_flank = ge_id.split('__')[-1]
        chromosome = ge_flank.split(':')[0]
        start, end = list(map(int, ge_flank.split(':')[-1].split('-')))
        length = end - start + 1
        nucleotide_sequence = hg38(chromosome, start - 1, size=length + 2)  # account for trinucleotide
        trinucleotides = [nucleotide_sequence[i:i + 3] for i in range(0, length)]

        # Compute bin trinucleotide probabilities from the reference sequence
        trin_probs = [trinuc_probs_dict.get(tn, 0) for tn in trinucleotides]
        prob_factor = 1 / sum(trin_probs)
        trin_probs_norm = [prob_factor * p for p in trin_probs]

        # Simulate mutations with trinucleotide sequence probability
        simulations = np.random.choice(range(0, length), size=(nsim, ge_mutations), p=trin_probs_norm, replace=True)

        # Save mutations per position on each iteration sepparately
        for i, iteration in enumerate(simulations):
            for position in iteration:
                total_simulations[i][position] += 1

    # Calculate fold change per simulation
    # Save simulated distribution of mutations per position
    expected = dict(zip(list(range(0, window_size)), [0] * window_size))
    expected_fc_list = []
    for iteration, hits_per_position in total_simulations.items():
        expected = Counter(expected) + Counter(hits_per_position)    # update expected
        unseen_positions = set(range(0, window_size)).difference(set(hits_per_position.keys()))
        unseen_hits = dict(zip(list(unseen_positions), [0] * len(unseen_positions)))
        hits_per_position.update(unseen_hits)
        _, muts_per_position = zip(*sorted(list(hits_per_position.items()), key=lambda x: x[0]))
        # Compute fold change of the simulation
        expected_fc_list += [compute_fold_change(muts_per_position, feature_size, flank_size, total_features)]
    # Expected mutations per position across simulations
    _, expected_muts_per_pos = zip(*sorted(list(expected.items()), key=lambda x: x[0]))

    # Test observed differs from expected
    obj = AnalyticalPvalue()
    obj.calculate_bandwidth(expected_fc_list)
    analytical_pvalue = obj.calculate(observed_fc)

    # Save results
    results = dict()
    results['observed_fc'] = observed_fc
    results['expected_fc'] = expected_fc_list
    results['pval'] = analytical_pvalue
    results['expected_muts'] = list(expected_muts_per_pos)
    with open(output_f, 'w') as ofd:
        json.dump(results, ofd)


if __name__ == '__main__':
    main()
