# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import scipy.spatial.distance as spdist

from skbio.util import get_rng


def _check_dist_metric(metric):
    """Validate distance metric."""
    half = False
    if isinstance(metric, str):
        if metric == "unitcorr":
            metric, half = spdist.correlation, True
        else:
            metric = getattr(spdist, metric)
    elif not callable(metric):
        raise ValueError("`metric` must be a string or callable.")
    return metric, half


def _check_shuffler(shuffler):
    """Validate sample shuffler."""
    if not callable(shuffler):
        shuffler = get_rng(shuffler).shuffle
    return shuffler
