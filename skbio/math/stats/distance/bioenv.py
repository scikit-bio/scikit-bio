# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.stats import spearmanr
from sklearn.preprocessing import scale


def bioenv(distance_matrix, data_frame, columns=None):
    """Find subset of variables maximally correlated with community distances.

    TODO: fill in description

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        TODO
    data_frame : pandas.DataFrame
        TODO
    columns : iterable of strs, optional
        TODO

    Returns
    -------
    pandas.DataFrame
        TODO

    Notes
    -----
    See [1]_ for the original BIO-ENV reference. The general algorithm and
    interface are similar to ``vegan::bioenv``, available in R's vegan package
    [2]_.

    References
    ----------
    .. [1] Clarke, K. R & Ainsworth, M. 1993. "A method of linking multivariate
       community structure to environmental variables". Marine Ecology Progress
       Series, 92, 205-219.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    """
    if columns is None:
        columns = data_frame.columns.values.tolist()

    if len(set(columns)) != len(columns):
        raise ValueError("Duplicate column names are not supported.")

    # TODO check for number of columns >= 1

    # TODO test for missing columns and/or ids
    # also add unit test to ensure differences in order don't affect results
    vars_df = data_frame.loc[distance_matrix.ids, columns]

    # TODO check for non-numeric columns

    # From http://stackoverflow.com/a/18017059
    # TODO make sure this doesn't modify the original df
    vars_df = pd.DataFrame(scale(vars_df), index=vars_df.index,
                           columns=vars_df.columns)
    dm_flat = distance_matrix.condensed_form()

    num_vars = len(columns)
    var_idxs = np.arange(num_vars)
    max_rhos = []
    max_vars = []
    for subset_size in range(1, num_vars + 1):
        # TODO performance test this way vs. generating all rhos in np array
        # and choosing max
        max_rho = (None, None)
        for subset_idxs in combinations(var_idxs, subset_size):
            vars_array = vars_df.iloc[:, subset_idxs].values
            vars_dm_flat = pdist(vars_array, metric='euclidean')
            rho = spearmanr(dm_flat, vars_dm_flat)[0]

            if max_rho == (None, None) or rho > max_rho[0]:
                max_rho = (rho, subset_idxs)
        max_rhos.append((subset_size, max_rho[0]))
        max_vars.append(', '.join([columns[i] for i in max_rho[1]]))

    return pd.DataFrame.from_records(max_rhos, index=max_vars,
                                     columns=('size', 'correlation'))
