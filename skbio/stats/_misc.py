# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import warnings

import numpy as np


def p_value_to_str(p_value, permutations):
    """Format p-value as a string with the correct number of decimals.

    .. note:: Deprecated in scikit-bio 0.2.1-dev
       ``p_value_to_str`` will be removed in scikit-bio 0.3.0.
       Permutation-based p-values in scikit-bio are calculated as
       ``(num_extreme + 1) / (num_permutations + 1)``, so it is impossible to
       obtain a p-value of zero. This function historically existed for
       correcting the number of digits displayed when obtaining a p-value of
       zero. Since this is no longer possible, this functionality will be
       removed.

    Number of decimals is determined by the number of permutations.

    Parameters
    ----------
    p_value : float or None
        p-value to convert to string.
    permutations : int
        Number of permutations used to calculate `p_value`.

    Returns
    -------
    str
        `p_value` formatted as a string with the correct number of decimals. If
        `p_value` is ``None`` or ``np.nan``, ``'N/A'`` is returned. If
        `permutations` is less than 10, a message stating insufficient number
        of permutations is returned.
    """
    warnings.warn(
        "skbio.stats.p_value_to_str is deprecated and will be removed in "
        "scikit-bio 0.3.0. There are no plans to provide a replacement for "
        "this functionality.", UserWarning)

    if p_value is None or np.isnan(p_value):
        result = 'N/A'
    elif permutations < 10:
        result = ('Too few permutations to compute p-value (permutations '
                  '= %d)' % permutations)
    else:
        decimal_places = int(np.log10(permutations + 1))
        result = ('%1.' + '%df' % decimal_places) % p_value

    return result
