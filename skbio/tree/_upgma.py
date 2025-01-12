# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from scipy.cluster.hierarchy import linkage
from skbio.tree import TreeNode
from ._utils import _check_dm


def upgma(dm, weighted=False):
    r"""Perform unweighted pair group method with arithmetic mean (UPGMA) or its
    weighted variant (WPGMA) for phylogenetic reconstruction.

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        The input distance matrix.
    weighted : bool, optional
        If True, WPGMA is performed instead of UPGMA. WPGMA is a variant of UPGMA
        which is unbiased towards the size of subtrees computed.

    Returns
    -------
    TreeNode
        A TreeNode object with estimated edge values.

    See Also
    --------
    nj

    Notes
    -----
    UPGMA (unweighted pair group method with arithmetic mean) is a simple hierarchical
    clustering method that iteratively groups proximal taxa or taxon groups to form a
    tree structure. A weighted variant is known as WPGMA, and both variants are due to
    Sokal and Michener [1]_.

    This function wraps SciPy's :func:`~scipy.cluster.hierarchy.linkage` function, with
    the ``method`` parameter set as "average" (UPGMA) or "weighted" (WPGMA). It takes a
    scikit-bio DistanceMatrix object and returns a scikit-bio TreeNode object.

    UPGMA creates a rooted and ultrametric tree -- all tips will have the same height
    (distance from the root node).

    References
    ----------
    .. [1] Sokal, R.R., & Michener, C.D. (1958). A statistical method for
       evaluating systematic relationships. University of Kansas science
       bulletin, 38, 1409-1438.

    Examples
    --------
    Define a distance matrix object for the taxa a, b, and c.

    >>> from skbio import DistanceMatrix

    >>> data = [[0, 1, 2],
    ...         [1, 0, 3],
    ...         [2, 3, 0]]
    >>> ids = list('abc')
    >>> dm = DistanceMatrix(data, ids)

    Construct a tree using UPGMA.

    >>> tree = upgma(dm)
    >>> print(tree.ascii_art())
              /-c
    ---------|
             |          /-a
              \--------|
                        \-b

    The tree also has estimated edge values assigned to each edge.

    >>> print(tree)
    (c:1.25,(a:0.5,b:0.5):0.75);
    <BLANKLINE>
    """
    _check_dm(dm)

    # If weighted is set to 'False', UPGMA is performed
    if weighted is False:
        linkage_matrix = linkage(
            dm.condensed_form(), method="average", metric="euclidean"
        )
    # Otherwise, WPGMA is performed
    else:
        linkage_matrix = linkage(
            dm.condensed_form(), method="weighted", metric="euclidean"
        )

    # Construct the TreeNode from the linkage matrix
    tree = TreeNode.from_linkage_matrix(linkage_matrix, dm.ids)

    return tree
