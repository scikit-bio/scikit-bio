# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, Union, Any, Callable, TYPE_CHECKING
from inspect import isfunction, getmodule

import numpy as np

from skbio.sequence import RNA
from skbio.stats.distance import DistanceMatrix
from skbio.sequence.distance import _check_seqtype, _char_hash, _char_freqs
import skbio.sequence.distance as sk_seqdist

if TYPE_CHECKING:  # pragma: no cover
    from skbio.alignment import TabularMSA


def align_dists(
    alignment: "TabularMSA",
    metric: Union[str, Callable],
    shared_by_all: Optional[bool] = True,
    **kwargs: Any,
) -> "DistanceMatrix":
    r"""Create a distance matrix from a multiple sequence alignment.

    .. versionadded:: 0.7.2

    This function calculates the distance between each pair of sequences based on the
    aligned sites, using a pre-defined distance metric under
    :mod:`skbio.sequence.distance` or a custom function. The resulting distance matrix
    can be used for phylogenetic tree reconstruction (see :mod:`skbio.tree`) or other
    types of distance-based analyses.

    Parameters
    ----------
    alignment : TabularMSA of shape (n_sequences, n_positions)
        Multiple sequence alignment.
    metric : str or callable
        The distance metric to apply to each pair of aligned sequences. See
        :mod:`skbio.sequence.distance` for available metrics. Or supply a custom
        function that takes two sequence objects and returns a number.
    shared_by_all : bool, optional
        Calculate the distance between each pair of sequences based on sites shared
        across all sequences (True, default), or shared between the current pair of
        sequences (False).
    kwargs : dict, optional
        Metric-specific parameters. Refer to the documentation of the chosen metric.

    Returns
    -------
    DistanceMatrix of shape (n_sequences, n_sequences)
        Distance matrix of aligned sequences.

    See Also
    --------
    skbio.sequence.distance
    skbio.stats.distance.DistanceMatrix

    Notes
    -----
    This function utilizes the preset metrics under :mod:`skbio.sequence.distance` to
    calculate pairwise sequence distances. These metrics have been optimized such that
    calling them through this function is significantly more efficient than directly
    calling them on every pair of sequences in an alignment. Custom functions may be
    provided but they may not have such acceleration.

    A "site" refers to a position in a sequence that has a valid character from an
    alphabet pre-defined by each metric. Typically, gaps are not considered as sites,
    and will be excluded from calculation. Some metrics are further restricted to
    certain characters, such as the four nucleobases or the 20 basic amino acids, while
    excluding degenerate and/or non-canonical characters.

    Sequences are filtered to aligned sites prior to calculation. This has two working
    modes: By default, alignment positions with at least one invalid character in any
    of the sequences are excluded from all sequences, leaving fully-aligned positions
    across the entire alignment. When ``shared_by_all`` is set to False, this filtering
    is applied to each pair of sequences, and positions with both characters valid are
    retained for calculation, even if there are invalid characters at the same position
    in other sequences.

    Some distance metrics require observed character frequencies (see the ``freqs``
    parameter of individual metric functions). They will be sampled from the entire
    alignment using all sites prior to filtering. Therefore, distances calculated using
    this function may be unequal to distances calculated by calling the metric function
    on each pair of sequences, which would sample character frequencies from that pair
    of sequences specifically.

    If you wish to apply a custom function on each pair of sequences, without filtering
    sites or considering global character frequencies, use
    ``DistanceMatrix.from_iterable(alignment, function)`` instead.

    Examples
    --------
    >>> from skbio.sequence import DNA
    >>> from skbio.alignment import TabularMSA, align_dists
    >>> msa = TabularMSA([
    ...     DNA('ATC-GTATCGG'),
    ...     DNA('ATGCG--CCGC'),
    ...     DNA('GTGCGTACGC-'),
    ... ], index=list("abc"))
    >>> dm = align_dists(msa, 'jc69')
    >>> print(dm)
    3x3 distance matrix
    IDs:
    'a', 'b', 'c'
    Data:
    [[ 0.          0.35967981  2.28339183]
     [ 0.35967981  0.          0.6354734 ]
     [ 2.28339183  0.6354734   0.        ]]

    """
    ids = alignment.index.astype(str)
    nseq, npos = alignment.shape

    # This is unusual. However, TabularMSA currently allows no sequence.
    if nseq == 0:
        raise ValueError("Alignment contains no sequence.")
    size = nseq * (nseq - 1) // 2

    # Also, TabularMSA allows zero-length sequences. Just return a matrix of NaN.
    if npos == 0:
        return DistanceMatrix(np.full(size, np.nan), ids=ids)

    # determine preset or custom function
    func, preset = _get_preset(metric)
    name = func.__name__

    # validate sequence type
    dtype = alignment.dtype
    seqtype = getattr(func, "_seqtype", None)
    _check_seqtype(name, dtype, seqtype)

    # get valid characters
    alphabet = getattr(func, "_alphabet", None)
    if alphabet is None and preset is False:
        alphabet = "nongap"
    valid = _char_hash(alphabet, dtype)

    # Create a 2D matrix of ASCII codes of all sequences.
    # TODO: This can be omitted after optimizing TabularMSA.
    seqs = np.vstack([seq._bytes for seq in alignment])

    # Create a 2D Boolean mask of valid sites.
    site_mat = None if alphabet is None else valid[seqs]

    # Get character frequencies from the entire alignment.
    if getattr(func, "_has_freqs", None) is True and "freqs" not in kwargs:
        kwargs["freqs"] = _char_freqs(seqs, valid)

    # Use a preset distance metric (efficient).
    if preset:
        func = getattr(sk_seqdist, "_" + name)

        # Delete positions that are not shared by all sequences (complete deletion). If
        # all positions are gone, return a NaN matrix.
        if shared_by_all and site_mat is not None:
            site_vec = np.all(site_mat, axis=0)
            seqs = seqs[:, site_vec]
            if seqs.shape[1] == 0:
                return DistanceMatrix(np.full(size, np.nan), ids=ids)
            site_mat = None

        # Otherwise, pass mask to the metric to trigger pairwise deletion.
        dm = func(seqs, site_mat, dtype, **kwargs)

    # Call a custom function on each pair of sequences (slow).
    else:
        dm = np.empty(size)

        # complete deletion
        if shared_by_all:
            site_vec = np.all(site_mat, axis=0)
            seqs = [seq[site_vec] for seq in alignment]
            pos = -1
            for i in range(nseq - 1):
                seq_i = seqs[i]
                for j in range(i + 1, nseq):
                    dm[pos + j] = func(seq_i, seqs[j], **kwargs)
                pos += nseq - i - 2

        # pairwise deletion
        else:
            pos = -1
            for i in range(nseq - 1):
                seq_i = alignment[i]
                sites_i = site_mat[i]
                for j in range(i + 1, nseq):
                    site_vec = sites_i & site_mat[j]
                    seq1 = seq_i[site_vec]
                    seq2 = alignment[j][site_vec]
                    dm[pos + j] = func(seq1, seq2, **kwargs)
                pos += nseq - i - 2

    return DistanceMatrix(dm, ids=ids)


def _get_preset(metric):
    """Validate a present distance metric and return the underlying function.

    Parameters
    ----------
    metric : str or callable
        The distance metric specified.

    Returns
    -------
    func : callable
        Function that calculates this metric.
    preset : bool
        Whether the metric is a preset one under `skbio.sequence.distance`.

    """
    if isinstance(metric, str):
        func = getattr(sk_seqdist, metric, None)
        if func is None or not isfunction(func) or not hasattr(func, "_is_metric"):
            raise ValueError(
                f"{metric!r} is not an available sequence distance metric name. Refer "
                "to `skbio.sequence.distance` for a list of available metrics."
            )
        preset = True

    elif callable(metric):
        func = metric
        if getmodule(func) is sk_seqdist:
            if not isfunction(func) or not hasattr(func, "_is_metric"):
                raise ValueError(
                    f"`{func.__name__}` is not a sequence distance metric. Refer to "
                    "`skbio.sequence.distance` for a list of available metrics."
                )
            preset = True
        else:
            preset = False

    else:
        raise TypeError("`metric` must be a function or a string.")

    return func, preset
