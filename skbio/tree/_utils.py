# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from ._exception import (
    NoLengthError,
    DuplicateNodeError,
    MissingNodeError,
    TreeError,
)


def _validate_dm(dm):
    if dm.shape[0] < 3:
        raise ValueError("Distance matrix must be at least 3x3 to generate a tree.")


def _validate_dm_and_tree(dm, tree):
    _validate_dm(dm)
    if frozenset(dm.ids) != tree.subset():
        raise ValueError("Inconsistent taxa between tree and distance matrix.")


def _validate_taxa_and_tree(taxa, tree, unique=True, rooted=False, lengths=False):
    """Validate taxa and tree prior to phylogenetic analyses.

    Parameters
    ----------
    taxa : sequence of str
        Input taxa.
    tree : TreeNode
        Input tree.
    unique : bool, optional
        If True, check if all taxa and tip names are unique.
    rooted : bool, optional
        If True, check if the tree is rooted.
    lengths : bool, optional
        If True, check if all non-root nodes in the tree have a branch length.

    Raises
    ------
    ValueError
        If there are duplicate taxa.
    TreeError
        If the tree is not rooted.
    NoLengthError
        If there are non-root nodes without a branch length.
    DuplicateNodeError
        If there are duplicate tip names.
    MissingNodeError
        If some taxa are not present as tip names in the tree.

    """
    # This function was moved from skbio.diversity._util and modified.

    taxon_set = set(taxa)
    if unique and len(taxa) != len(taxon_set):
        raise ValueError("All taxa must be unique.")

    # The following code checks if the tree is rooted.
    # There was a comment in the original code: "this is an imperfect check for "
    # "whether the tree is rooted or not. can this be improved?"
    # This check could be simplified as `not tree._is_rooted()`. It is kept as the
    # original form for compatibility.
    if rooted and len(tree.root().children) > 2:
        raise ValueError("The tree must be rooted.")

    tip_names = []
    tip_names_append = tip_names.append
    for node in tree.postorder(include_self=False):
        if lengths and node.length is None:
            raise ValueError(  # NoLengthError
                "All non-root nodes in the tree must have a branch length."
            )
        if not node.children:
            tip_names_append(node.name)

    tip_name_set = set(tip_names)
    if unique and len(tip_names) != len(tip_name_set):
        raise DuplicateNodeError("All tip names in the tree must be unique.")

    if missing := taxon_set - tip_name_set:
        raise MissingNodeError(
            f"{len(missing)} taxa are not present as tip names in the tree."
        )


def _tree_to_lnkmat(tree, taxa=None, out=None):
    """Convert a bifurcating tree into a linkage matrix.

    This function creates a matrix that resembles SciPy's linkage matrix, but only
    retains the first two columns (left and right children per internal node). This is
    because the two children may have uneven branch lengths, which is unlike SciPy's
    `linkage` function's output.

    Parameters
    ----------
    tree : TreeNode
        Input tree.
    taxa : sequence of str, optional
        Taxon name in order. If provided, the output will adopt this order. Otherwise,
        will use tip names in postorder.
    out : ndarray of shape (n - 1, 2), optional
        Location to write output linkage matrix into. _n_ is the number of taxa. If not
        provided, an array of `np.intp` type will be created.

    Returns
    -------
    ndarray of shape (n - 1, 2)
        Output linkage matrix.

    Notes
    -----
    Order of (left, right) children of each internal node is preserved.

    """
    # TODO: Parse uneven branch lengths, as `_nj._tree_from_linkmat` does.
    if taxa is None:
        taxa = [tip.name for tip in tree.tips()]
    n = len(taxa)

    index = dict(zip(taxa, range(n)))  # taxon-to-index mapping
    if len(index) < n:
        raise ValueError("Taxa contain duplicates.")

    if out is None:
        out = np.empty((n - 1, 2), dtype=np.intp)
    elif out.shape != (n - 1, 2):
        raise ValueError(f"`out` must have a shape of ({n - 1}, 2).")

    nodes = {}  # node-to-index mapping
    nodes_pop = nodes.pop

    i = 0
    for node in tree.postorder(include_self=True):
        if not node.children:
            try:
                nodes[node] = index[node.name]
            except KeyError:
                raise ValueError(f"Tip name {node.name!r} is absent from taxa.")
        else:
            try:
                a, b = node.children
            except ValueError:
                raise ValueError("Tree must be strictly bifurcating.")
            try:
                out[i, 0] = nodes_pop(a)
                out[i, 1] = nodes_pop(b)
            except IndexError:  # when taxa are provided but tree has duplicates
                raise ValueError("Tree contains duplicate tip names.")
            nodes[node] = n + i
            i += 1

    if i != n - 1:
        raise ValueError("One or more taxa are absent from the tree.")
    if len(np.unique(out)) < out.size:
        raise ValueError("Tree contains duplicate tip names.")

    return out
