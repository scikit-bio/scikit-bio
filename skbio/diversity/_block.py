# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.diversity._driver import partial_beta_diversity
from skbio.stats.distance import DistanceMatrix
from skbio.diversity._util import _validate_counts_matrix


def _generate_id_blocks(ids, k=64):
    """Generate blocks of IDs that map into a DistanceMatrix.

    Parameters
    ----------
    ids : Iterable object
        An iterable of IDs of whatever type.
    k : int, optional
        The size of a block to generate IDs for, defaults to 64.

    Notes
    -----
    This method is intended to facilitate partial beta diversity calculations.
    Blocks of IDs are generated from the upper triangle of the subsequent
    distance matrix. For instance, given the following distance matrix with
    IDs {A, B, C, D, E}:

      A B C D E
    A 0 # # # #
    B # 0 # # #
    C # # 0 # #
    D # # # 0 #
    E # # # # 0

    The goal of this method is to generate tuples of IDs of at most size k over
    the upper triangle which correspond to blocks of the matrix to compute. IDs
    are remapped as well into integers to facilitate downstream indexing.
    Given k=3, the following ID tuples would be generated:

    ((0, 1, 2), (0, 1, 2))
    ((0, 1, 2), (3, 4))
    ((3, 4), (3, 4))

    This method is not responsible for describing which specific pairs of IDs
    are to be computed, only the subset of the matrix of interest.

    Returns
    -------
    tuple of 1D np.array
        Index 0 contains the row IDs, and index 1 contains the column IDs

    """
    n = len(ids)
    ids_idx = np.arange(n)

    for row_start in range(0, n, k):
        for col_start in range(row_start, n, k):
            row_ids = ids_idx[row_start : row_start + k]
            col_ids = ids_idx[col_start : col_start + k]

            yield (row_ids, col_ids)


def _block_party(counts=None, row_ids=None, col_ids=None, **kwargs):
    """Subset counts to relevant rows and columns.

    Parameters
    ----------
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of taxa in a given sample.
    row_ids : 1D np.ndarray of int
        Block row IDs to keep in the counts matrix.
    col_ids : 1D np.ndarray of int
        Block column IDs to keep in the counts matrix. Note, these correspond
        to rows in the counts matrix, but columns in a subsequent distance
        matrix.
    kwargs : dict
        Keyword arguments containing information about the block to compute.

    Returns
    -------
    dict
        kwargs that describe the block to compute. A filtered ``counts`` matrix
        is stored in kwargs. If applicable, a filtered ``tree`` and ``taxa``
        are also stored.

    """
    ids_to_keep = np.unique(np.hstack([row_ids, col_ids]))

    # create a view of the relevant samples
    counts_block = counts[ids_to_keep]

    # remove from the block any empty observations
    # NOTE: this will perform an implicit copy
    nonzero_cols = (counts_block != 0).any(axis=0)
    counts_block = counts_block[:, nonzero_cols]

    kwargs["counts"] = counts_block
    kwargs["ids"] = ids_to_keep

    if "tree" in kwargs and "taxa" in kwargs:
        kwargs["taxa"] = np.asarray(kwargs["taxa"])[nonzero_cols]
        kwargs["tree"] = kwargs["tree"].shear(kwargs["taxa"])

    return kwargs


def _pairs_to_compute(rids, cids):
    """Determine the pairs of samples to compute distances between.

    Parameters
    ----------
    rids : Iterable
        The row IDs in the partial pairwise computation.
    cids : Iterable
        The column IDs in the partial pairwise computation.

    Raises
    ------
    ValueError
        When determining ID pairs for blocks that fall outside of the diagonal
        of the resulting distance matrix, if a pair corresponds to the lower
        triangle, complain loudly.

    Returns
    -------
    list of tuple
        The ID pairs to compute distances between.

    """
    # if identical, gather the upper triangle
    if len(rids) == len(cids) and (rids == cids).all():
        return [(i, j) for idx, i in enumerate(rids) for j in rids[idx + 1 :]]

    # otherwise, grab pairwise combinations disregarding the diagonal
    else:
        if set(rids).intersection(set(cids)):
            raise ValueError("Attempting to compute a lower triangle")
        return [(i, j) for i in rids for j in cids if i != j]


def _block_kwargs(**kwargs):
    """Construct arguments describing a block to compute.

    Returns
    -------
    dict
        The parameters for the block of the distance matrix to compute.

    """
    valid_block_keys = {
        "counts",
        "ids",
        "tree",
        "taxa",
        "metric",
        "id_pairs",
        "validate",
        "otu_ids",
    }
    for row_ids, col_ids in _generate_id_blocks(kwargs["ids"], kwargs["k"]):
        id_pairs = _pairs_to_compute(row_ids, col_ids)
        if id_pairs:
            kw = {k: v for k, v in kwargs.items() if k in valid_block_keys}
            kw["id_pairs"] = id_pairs
            kw["row_ids"] = row_ids
            kw["col_ids"] = col_ids

            yield kw


def _block_compute(**kwargs):
    """Compute a block within the resulting distance matrix.

    Notes
    -----
    This method encapsulates the two expensive operations to perform for each
    block, namely, the "shearing" of the phylogenetic tree to correspond to
    only the taxa of interest, and the actual beta diversity calculations.

    Returns
    -------
    DistanceMatrix

    """
    block_kw = _block_party(**kwargs)

    return partial_beta_diversity(**block_kw)


def _map(func, kw_gen):
    """Map a function over arguments.

    Notes
    -----
    builtin map does not allow for mapping with kwargs.

    Parallel uses of block decomposition will likely replace this method with
    one which can distribute compute.

    """
    for kwargs in kw_gen:
        yield func(**kwargs)


def _reduce(blocks):
    """Reduce an iterable of partial distance matrices into a full matrix.

    Note, the reduce doesn't actually care about what pairs are computed
    so if a distance between pairs exists multiple times, it'll get
    added. as such, this reduction is only safe to perform if by
    the block_beta_diversity method which assures that distances are not
    computed multiple times.
    """
    all_blocks = list(blocks)

    # Determine the maximum integer ID observed in the blocks. There exists a
    # 1-1 mapping between the integer ID and a sample ID. We increment by 1
    # as the integer ID space begins with zero, and we'll be using this value
    # to determine the size of the resulting full distance matrix.
    n_ids = max(map(lambda x: max(x.ids), all_blocks)) + 1

    mat = np.zeros((n_ids, n_ids), dtype=float)

    # TODO: something smarter.
    for block in all_blocks:
        n_blk_ids = len(block.ids)

        # get the corresponding coordinates in the master matrix
        master_idx = [
            (i, j) for row, i in enumerate(block.ids) for j in block.ids[row + 1 :]
        ]

        # get the corresponding coordinates within the current block
        block_idx = [
            (i, j)
            for row, i in enumerate(range(n_blk_ids))
            for j in range(row + 1, n_blk_ids)
        ]

        for (m_i, m_j), (b_i, b_j) in zip(master_idx, block_idx):
            mat[m_i, m_j] += block.data[b_i, b_j]

    return DistanceMatrix(mat + mat.T, list(range(n_ids)))


def block_beta_diversity(
    metric, counts, ids, validate=True, k=64, reduce_f=None, map_f=None, **kwargs
):
    """Perform a block-decomposition beta diversity calculation.

    Parameters
    ----------
    metric : str or callable
        The pairwise distance function to apply. If ``metric`` is a string, it
        must be resolvable by scikit-bio (e.g., UniFrac methods), or must be
        callable.
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of taxa in a given sample.
    ids : iterable of strs
        Identifiers for each sample in ``counts``.
    validate : bool, optional
        See ``skbio.diversity.beta_diversity`` for details.
    reduce_f : function, optional
        A method to reduce `PartialDistanceMatrix` objects into a single
        `DistanceMatrix`. The expected signature is:

            `f(Iterable of DistanceMatrix) -> DistanceMatrix`

        Note, this is the reduce within a map/reduce.
    map_f: function, optional
        A method that accepts a `_block_compute`. The expected signature is:

            `f(**kwargs) -> DistanceMatrix`

        NOTE: ipyparallel's `map_async` will not work here as we need to be
        able to pass around `**kwargs``.
    k : int, optional
        The blocksize used when computing distances
    kwargs : kwargs, optional
        Metric-specific parameters.

    Returns
    -------
    DistanceMatrix
        A distance matrix relating all samples represented by counts to each
        other.

    Notes
    -----
    This method is designed to facilitate computing beta diversity in parallel.
    In general, if you are processing a few hundred samples or less, then it is
    likely the case that `skbio.diversity.beta_diversity` will be faster. The
    original need which motivated the development of this method was processing
    the Earth Microbiome Project [1]_ dataset which at the time spanned over
    25,000 samples and 7.5 million open reference taxa.

    See Also
    --------
    skbio.diversity.beta_diversity
    skbio.diversity.partial_beta_diversity

    References
    ----------
    .. [1] http://www.earthmicrobiome.org/

    """
    if validate:
        counts = _validate_counts_matrix(counts, ids=ids)

    if reduce_f is None:
        reduce_f = _reduce

    if map_f is None:
        map_f = _map

    # The block method uses numeric IDs to take advantage of fancy indexing
    # with numpy.
    tmp_ids = np.arange(len(counts))
    kwargs["ids"] = tmp_ids

    kwargs["metric"] = metric
    kwargs["counts"] = counts
    kwargs["k"] = k
    kwargs["validate"] = False  # we've already validated if necessary

    dm = reduce_f(map_f(_block_compute, _block_kwargs(**kwargs)))
    dm.ids = ids

    return dm
