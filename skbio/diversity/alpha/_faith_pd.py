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
from skbio.diversity._validate import (_validate_counts_vector,
                                       _validate_otu_ids_and_tree,
                                       _vectorize_counts_and_tree)


def _faith_pd(counts_by_node, branch_lengths):
    observed_nodes = np.where(counts_by_node > 0, 1, 0)
    return (branch_lengths * observed_nodes).sum()


@experimental(as_of="0.4.0-dev")
def faith_pd(counts, otu_ids, tree, validate=True):
    """ Compute Faith's phylogenetic diversity metric (PD)

    Parameters
    ----------
    counts : 1-D array_like, int
        Vectors of counts/abundances of OTUs for one sample.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``counts``.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.
    validate: bool, optional
        If `False`, validation of the input won't be performed. This step can
        be slow, so if validation is run elsewhere it can be disabled here.
        However, invalid input data can lead to invalid results or error
        messages that are hard to interpret, so this step should not be
        bypassed if you're not certain that your input data are valid. See
        Notes for the description of what validation entails so you can
        determine if you can safely disable validation.

    Returns
    -------
    float
        The phylogenetic diversity (PD) of the samples.

    Raises
    ------
    ValueError, MissingNodeError, DuplicateNodeError
        If validation fails (see description of validation in Notes). Exact
        error will depend on what was invalid.

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

    Validation of input data confirms the following:
     * ``counts`` data can be safely cast to integers
     * there are no negative values in ``counts``
     * ``counts`` has the correct number of dimensions
     * ``otu_ids`` does not contain duplicate values
     * ``len(counts)`` is equal to ``len(otu_ids)``
     * ``tree`` is rooted
     * ``tree`` has more than one node
     * all nodes in ``tree`` except for the root node have branch lengths
     * all tip names in ``tree`` are unique
     * all ``otu_ids`` correspond to tip names in ``tree``

    References
    ----------
    .. [1] Faith, D. P. Conservation evaluation and phylogenetic diversity.
       Biol. Conserv. (1992).

    """
    counts_by_node, branch_lengths = _setup_faith_pd(
        counts, otu_ids, tree, validate, single_sample=True)

    return _faith_pd(counts_by_node, branch_lengths)


def _setup_faith_pd(counts, otu_ids, tree, validate, single_sample):
    if validate:
        if single_sample:
            # only validate count if operating in single sample mode, they
            # will have already been validated otherwise
            counts = _validate_counts_vector(counts)
            _validate_otu_ids_and_tree(counts, otu_ids, tree)
        else:
            _validate_otu_ids_and_tree(counts[0], otu_ids, tree)

    counts_by_node, tree_index, branch_lengths = \
        _vectorize_counts_and_tree(counts, otu_ids, tree)

    return counts_by_node, branch_lengths
