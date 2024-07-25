# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from scipy.cluster.hierarchy import average, weighted
from skbio import TreeNode
from skbio.stats.distance import DistanceMatrix


def upgma(dm: DistanceMatrix, weighted=False) -> TreeNode:
    r"""
    This function implements the Unweighted Pair Group Method with Arithmetic
    mean (UPGMA) algorithm as well as the weighted variant (WPGMA) for phylogenetic
    reconstruction.

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        The input distance matrix.
    weighted : bool, optional
        If 'True', WPGMA is performed instead of UPGMA. WPGMA is a variant of
        UPGMA which is unbiased towards the size of subtrees computed.

    Returns
    -------
    TreeNode
        A TreeNode object with estimated edge values

    Notes
    -----
    UPGMA is a hierarchical clustering method appearing as the 'average' function
    in the SciPy package, where the linkage matrix produced by 'average' is used
    to construct a TreeNode object. A weighted variant known as WPGMA is attributed
    to Sokal and Michener [1]_.

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

    """
    # Ensure the input is a DistanceMatrix object
    if not isinstance(dm, DistanceMatrix):
        raise ValueError("Input must be a DistanceMatrix object.")

    # If weighted is set to 'False', UPGMA is performed
    if weighted is False:
        linkage_matrix = average(dm.condensed_form())
    # Otherwise, WPGMA is performed
    else:
        linkage_matrix = weighted(dm.condensed_form())

    # Construct the TreeNode from the linkage matrix
    tree = TreeNode.from_linkage_matrix(linkage_matrix, dm.ids)

    return tree
