# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from scipy.cluster.hierarchy import linkage
from skbio.tree import TreeNode
from skbio.stats.distance import DistanceMatrix


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
    UPGMA is a hierarchical clustering method appearing as the `average` function
    in the SciPy package, where the linkage matrix produced by `average` is used
    to construct a TreeNode object. A weighted variant is known as WPGMA, and both
    variants are due to Sokal and Michener [1]_.

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
    # Ensure the input is a DistanceMatrix object
    if not isinstance(dm, DistanceMatrix):
        raise ValueError("Input must be a DistanceMatrix object.")

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
