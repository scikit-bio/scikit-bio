"""
Miscellaneous statistics utilities (:mod:`skbio.math.stats.misc`)
=================================================================

.. currentmodule:: skbio.math.stats.misc

This module contains miscellaneous functionality that is generally useful when
working with statistical methods or their results in scikit-bio.

Functions
---------

.. autosummary::
   :toctree: generated/

   p_value_to_str

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np


def p_value_to_str(p_value, permutations):
    """Format p-value as a string with the correct number of decimals.

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
    if p_value is None or np.isnan(p_value):
        result = 'N/A'
    elif permutations < 10:
        result = ('Too few permutations to compute p-value (permutations '
                  '= %d)' % permutations)
    else:
        decimal_places = int(np.log10(permutations + 1))
        result = ('%1.' + '%df' % decimal_places) % p_value

    return result
