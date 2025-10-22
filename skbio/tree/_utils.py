# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------


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
