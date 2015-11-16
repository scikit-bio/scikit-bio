# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

from skbio.util._decorator import experimental
from skbio.diversity._driver import (_validate_counts_vector,
                                     _validate_otu_ids_and_tree,
                                     _vectorize_counts_and_tree)


def _faith_pd(counts_by_node, branch_lengths):
    counts_by_node = np.where(counts_by_node > 0, 1, 0)
    result = (branch_lengths * counts_by_node).sum()
    return result


@experimental(as_of="0.4.0-dev")
def faith_pd(counts, otu_ids, tree, validate=True,):
    """ Compute Faith's phylogenetic diversity metric (PD)

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``counts``.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.
    validate: bool, optional
        If `False`, validation of the input won't be performed. This step can
        be slow, so if validation is run elsewhere it can be disabled here.
        However, invalid input data can lead to invalid results, so this step
        should not be bypassed all together.

    Returns
    -------
    float
        The phylogenetic diversity (PD) of the samples.

    Raises
    ------
    ValueError
        If ``counts`` and ``otu_ids`` are not equal in length.
    MissingNodeError
        If an OTU id is provided that does not correspond to a tip in the
        tree.

    Notes
    -----
    Faith's phylogenetic diversity, often referred to as PD, was originally
    described in [1]_.

    This implementation differs from that in PyCogent (and therefore QIIME
    versions less than 2.0.0) by imposing a few additional restrictions on the
    inputs. First, the input tree must be rooted. In PyCogent, if an unrooted
    tree was provided that had a single trifurcating node (a newick convention
    for unrooted trees) that node was considered the root of the tree. Next,
    all OTU IDs must be tips in the tree. PyCogent would silently ignore OTU
    IDs that were not present the tree. To reproduce Faith PD results from
    PyCogent with scikit-bio, ensure that your PyCogent Faith PD calculations
    are performed on a rooted tree and that all OTU IDs are present in the
    tree.

    References
    ----------
    .. [1] Faith, D. P. Conservation evaluation and phylogenetic diversity.
       Biol. Conserv. (1992).

    """
    if validate:
        counts = _validate_counts_vector(counts)
        _validate_otu_ids_and_tree(counts, otu_ids, tree)

    counts = np.asarray(counts)
    otu_ids = np.asarray(otu_ids)

    if counts.sum() == 0:
        return 0.0

    counts_by_node, tree_index, branch_lengths = \
        _vectorize_counts_and_tree(counts, otu_ids, tree)

    return _faith_pd(counts_by_node, branch_lengths)
