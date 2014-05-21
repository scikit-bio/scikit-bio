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
from scipy.stats import spearmanr

from skbio.core.distance import DistanceMatrix


def bioenv(distance_matrix, data_frame, columns):
    """BIO-ENV statistical method executor.

    TODO: fill in description

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
    cats = columns
    dm = distance_matrix
    dm_flat = dm.condensed_form()

    row_count = dm.shape[0]
    col_count = len(cats)
    sum = 0
    stats = [(-777777777, '') for c in range(col_count + 1)]
    for i in range(1, col_count + 1):
        combo = combinations([j for j in range(1, col_count + 1)], i)

        for element in combo:
            cat_mat = _make_cat_mat(cats, element, dm, data_frame)
            cat_dm = _derive_euclidean_dm(cat_mat, row_count, dm.ids)
            cat_dm_flat = cat_dm.condensed_form()
            r = spearmanr(dm_flat, cat_dm_flat)[0]
            if r > stats[i - 1][0]:
                stats[i - 1] = (r, ','.join(str(s) for s in element))

    res = {}
    res['method_name'] = 'BEST'
    res['num_vars'] = col_count
    res['vars'] = ['%s = %d' % (name, val + 1)
                   for val, name in enumerate(cats)]
    res['rho_vals'] = stats[:-1]

    return res

def _derive_euclidean_dm(cat_mat, dim, ids):
    """Returns an n x n, euclidean distance matrix, where n = len(cats)."""
    res_mat = []

    for i in range(dim):
        res_mat.append([0 for k in range(dim)])
        for j in range(i):
            res_mat[i][j] = _vector_dist(cat_mat[i], cat_mat[j])
            res_mat[j][i] = res_mat[i][j]

    return DistanceMatrix(res_mat, ids)

def _vector_dist(vec1, vec2):
    """Calculates the Euclidean distance between two vectors."""
    return np.sqrt(sum([(float(v1) - float(v2)) ** 2 for v1, v2 in
                     zip(vec1, vec2)]))

def _make_cat_mat(cats, combo, dm, df):
    """Returns a matrix with columns pulled from category values.

    Returns a matrix with len(ids) rows of columns pulled from
    category values, the number of columns for each category is
    determined by the current combination (combo).
    """
    res = []
    for i in combo:
        row = []
        for id_ in dm.ids:
            row.append(df[cats[i - 1]][id_])
        res.append(row)
    return zip(*res)
