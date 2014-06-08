# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np
from scipy.stats import pearsonr, spearmanr

from skbio.core.distance import DistanceMatrix

def mantel(x, y, method='pearson', permutations=999, alternative='twosided'):
    if method == 'pearson':
        corr_func = pearsonr
    elif method == 'spearman':
        corr_func = spearmanr
    else:
        raise ValueError("Invalid correlation method '%s'." % method)

    if permutations < 0:
        raise ValueError("Number of permutations must be greater than or "
                         "equal to zero.")
    if alternative not in ('twosided', 'greater', 'less'):
        raise ValueError("Invalid alternative hypothesis '%s'." % alternative)

    x = DistanceMatrix(x)
    y = DistanceMatrix(y)

    # TODO: test size >= 3
    if x.shape != y.shape:
        raise ValueError("Distance matrices must have the same shape.")

    x_flat = x.condensed_form()
    y_flat = y.condensed_form()

    orig_stat = corr_func(x_flat, y_flat)[0]

    if permutations == 0:
        p_value = np.nan
    else:
        size = x.shape[0]
        better = 0
        for i in range(permutations):
            perm = DistanceMatrix(_permute_2d(x, np.random.permutation(size)))
            perm_flat = perm.condensed_form()
            r = corr_func(perm_flat, y_flat)[0]

            if alternative == 'twosided':
                if abs(r) >= abs(orig_stat):
                    better += 1
            else:
                if ((alternative == 'greater' and r >= orig_stat) or
                    (alternative == 'less' and r <= orig_stat)):
                    better += 1
        p_value = (better + 1) / (permutations + 1)

    return orig_stat, p_value


def _permute_2d(x, p):
    """Performs 2D permutation of matrix x according to p."""
    return x[p][:, p]
