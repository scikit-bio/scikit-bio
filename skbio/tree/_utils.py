# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------


def _check_dm(dm):
    if dm.shape[0] < 3:
        raise ValueError("Distance matrix must be at least 3x3 to generate a tree.")


def _check_dm_tree(dm, tree):
    _check_dm(dm)
    if frozenset(dm.ids) != tree.subset():
        raise ValueError("Inconsistent taxa between tree and distance matrix.")
