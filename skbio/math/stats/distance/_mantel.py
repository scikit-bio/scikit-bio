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
    # TODO: test size >= 3

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

    x = np.asarray(x)
    y = np.asarray(y)

    if x.shape != y.shape:
        raise ValueError("Distance matrices must have the same shape.")
    # TODO assert square
    if not (_is_symmetric_and_hollow(x) and _is_symmetric_and_hollow(y)):
        raise ValueError("Distance matrices must be symmetric and hollow.")

    # Get a flattened list of lower-triangular matrix elements (excluding the
    # diagonal) in column-major order. Use these values to calculate the
    # correlation statistic.
    x_flat, y_flat = _flatten_lower_triangle(x), _flatten_lower_triangle(y)
    orig_stat = corr_func(x_flat, y_flat)[0]

    size = len(x)
    better = 0
    perm_stats = []
    for i in range(permutations):
        perm = _permute_2d(x, np.random.permutation(size))
        perm_flat = _flatten_lower_triangle(perm)
        r = corr_func(perm_flat, y_flat)[0]

        if alternative == 'twosided':
            if abs(r) >= abs(orig_stat):
                better += 1
        else:
            if ((alternative == 'greater' and r >= orig_stat) or
                (alternative == 'less' and r <= orig_stat)):
                better += 1
        perm_stats.append(r)

    return orig_stat, (better + 1) / (permutations + 1)


def _is_symmetric_and_hollow(x):
    return (x.T == x).all() and (np.trace(x) == 0)


def _flatten_lower_triangle(x):
    """Returns a list containing the flattened lower triangle of the matrix.

    The returned list will contain the elements in column-major order. The
    diagonal will be excluded.

    """
    x = np.asarray(x)
    flattened = []
    for col_num in range(x.shape[1]):
        for row_num in range(x.shape[0]):
            if col_num < row_num:
                    flattened.append(x[row_num][col_num])
    return flattened


def _permute_2d(x, p):
    """Performs 2D permutation of matrix x according to p."""
    return x[p][:, p]
