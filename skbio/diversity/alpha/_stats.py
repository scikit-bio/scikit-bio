# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function


def phylogenetic_diversity(u, otu_ids, tree):
    """ Compute Faith's phylogenetic diversity metric (PD)

    Parameters
    ----------
    u: list, np.array
        Vector of counts of OTUs.
    otu_ids: list, np.array
        Vector of OTU ids corresponding to tip names in ``tree``. Must be the
        same length as ``u``.
    tree: skbio.TreeNode
        Tree relating the OTUs in otu_ids. The set of tip names in the tree can
        be a superset of ``otu_ids``, but not a subset.

    Returns
    -------
    float
        The phylogenetic diversity (PD) of the samples.

    Raises
    ------
    ValueError
        If ``u`` and ``otu_ids`` are not equal in length.
    MissingNodeError
        If an OTU id is provided that does not correspond to a tip in the
        tree.

    Notes
    -----
    Phylogenetic diversity, often referred to as PD, was originally described
    in [1]_.

    References
    ----------
    .. [1] Faith, D. P. Conservation evaluation and phylogenetic diversity.
       Biol. Conserv. (1992).

    """
    observed_otus = {o: u for o, u in zip(otu_ids, u) if u >= 1}
    observed_nodes = tree.observed_node_counts(observed_otus)
    result = sum(o.length for o in observed_nodes if o.length is not None)
    return result
