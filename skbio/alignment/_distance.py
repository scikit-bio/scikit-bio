# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, Union, Iterable, Any, Callable, TYPE_CHECKING
from itertools import combinations
from inspect import isfunction

import numpy as np

from skbio.stats.distance import DistanceMatrix

if TYPE_CHECKING:  # pragma: no cover
    from skbio.alignment import TabularMSA


def align_dists(
    alignment: "TabularMSA",
    metric: Union[str, Callable],
    ids: Optional[Union[str, Iterable[Any]]] = None,
    shared_by_all: Optional[bool] = True,
    **kwargs: Any,
) -> "DistanceMatrix":
    r"""Create a distance matrix from a multiple sequence alignment.

    .. versionadded:: 0.7.1

    This function calculates the distance between each pair of sequences based on the
    aligned sites and a given distance metric. The resulting distance matrix can be
    used for phylogenetic tree reconstruction (see :mod:`skbio.tree`) or other types of
    distance-based analyses.

    Parameters
    ----------
    alignment : TabularMSA of shape (n_sequences, n_positions)
        Multiple sequence alignment.
    metric : str or callable
        The distance metric to apply to each pair of aligned sequences. See
        :mod:`skbio.sequence.distance` for available metrics. Passing metric as a
        string is preferable as this often results in an optimized version of the
        metric being used.
    ids : str or array_like of shape (n_sequences,), optional
        Unique identifiers of sequences in the alignment. Can be a list of strings, or
        a key to the ``metadata`` property of each sequence. If not provided, the
        metadata key ``id`` will be used if every sequence has it. Otherwise,
        incremental integer strings "1", "2", "3",... will be used.
    shared_by_all : bool, optional
        Calculate the distance between each pair of sequences based on sites shared
        across all sequences (True, default), or shared between the current pair of
        sequences (False). A "site" refers to a position in a sequence that has a
        definite character (i.e., not a gap or a degenerate character).
    kwargs : dict, optional
        Metric-specific parameters. Refer to the documentation of the chosen metric.

    Returns
    -------
    DistanceMatrix of shape (n_sequences, n_sequences)
        Distance matrix of aligned sequences.

    See Also
    --------
    skbio.sequence.distance

    Notes
    -----
    This function always strips positions with gap or degenerate characters from the
    aligned sequences. If you intend to compute distances using unstripped sequences,
    consider using ``DistanceMatrix.from_iterable(alignment, metric)``.

    """
    n_seqs, n_pos = alignment.shape

    # Identify aligned sites (definite characters).
    sites = [seq.definites() for seq in alignment]

    # Use a preset distance metric (efficient).
    if isinstance(metric, str):
        import skbio.sequence.distance

        func = getattr(skbio.sequence.distance, metric, None)
        if func is None or not isfunction(func) or not hasattr(func, "_is_metric"):
            raise ValueError(
                f'"{metric}" is not an available sequence distance metric name. '
                "Refer to `skbio.sequence.distance` for a list of available "
                "metrics."
            )
        func_name = "_" + metric
        if shared_by_all is False:
            func_name += "_pair"
        func = getattr(skbio.sequence.distance, func_name)

        # Filter sequences by a given alphabet.
        alphabet = func._alphabet

        if alphabet is not None:
            if alphabet in ("nongap", "definite", "canonical"):
                valid = getattr(alignment[0], f"_{alphabet}_hash")
            else:
                encoded = _encode_alphabet(alphabet)
                valid = np.zeros((Sequence._num_ascii_codes,), dtype=bool)
                valid[encoded] = True

            if equal:
                pos = valid[seq1._bytes] & valid[seq2._bytes]
                seq1, seq2 = seq1[pos], seq2[pos]
            else:
                seq1, seq2 = seq1[valid[seq1._bytes]], seq2[valid[seq2._bytes]]

        # Create 2D arrays of sequences (ASCII codes) and sites (Boolean mask).
        # This can be omitted after optimizing TabularMSA.
        sites = np.vstack(sites)
        seqs = [seq._bytes for seq in alignment]
        seqs = np.vstack(seqs)

        # complete deletion
        if shared_by_all:
            sites = np.all(sites, axis=0)
            seqs = seqs[:, sites]
            dm = func(seqs, **kwargs)

        # pairwise deletion
        else:
            dm = func(seqs, sites, **kwargs)

    # Call a custom function on each pair of sequences (slow).
    elif callable(metric):
        dm = np.empty((n_seqs * (n_seqs - 1) // 2,))

        # complete deletion
        if shared_by_all:
            sites = np.all(np.vstack(sites), axis=0)
            seqs = [seq[sites] for seq in alignment]
            pos = -1
            for i in range(n_seqs - 1):
                for j in range(i + 1, n_seqs):
                    dm[pos + j] = metric(seqs[i], seqs[j], **kwargs)
                pos += n_seqs - i - 2

        # pairwise deletion
        else:
            pos = -1
            for i in range(n_seqs - 1):
                for j in range(i + 1, n_seqs):
                    sites_ = sites[i] & sites[j]  ## will stripe twice
                    seq1 = alignment[i][sites_]
                    seq2 = alignment[j][sites_]
                    dm[pos + j] = metric(seq1, seq2, **kwargs)
                pos += n_seqs - i - 2
    else:
        raise TypeError("`metric` must be a function or a string.")

    return DistanceMatrix(dm, ids=alignment.index.astype(str))


def _valid_chars(alphabet, seqs):
    """Create a mask of valid characters for distance calculation.

    Parameters
    ----------
    alphabet : str, 1D array_like, {'nongap', 'definite', 'canonical'}
        An alphabet of valid characters to be considered by the metric.

    Returns
    -------
    ndarray of bool of shape (n_sequences, n_positions)
        Boolean mask of valid characters (True).

    See Also
    --------
    skbio.sequence.distance._metric_specs

    """
    from skbio.sequence._alphabet import _encode_alphabet

    n_seqs = len(seqs)
    seq1 = seqs[0]
    L = len(seq1)
    if alphabet in ("nongap", "definite", "canonical"):
        valid = getattr(seq1, f"_{alphabet}_hash")
    else:
        encoded = _encode_alphabet(alphabet)
        valid = np.zeros((seq1._num_ascii_codes,), dtype=bool)
        valid[encoded] = True

    mask = np.vstack([valid[seq._bytes] for seq in seqs])


def _ids_from_msa(msa, ids):
    """Get sequence IDs from an alignment."""
    key = None
    if ids is None:
        key = "id"
    elif isinstance(ids, str):
        key = ids

    if key is not None:
        try:
            ids = [seq.metadata[key] for seq in msa]
        except KeyError:
            if ids is not None:
                raise KeyError(f"Metadata key '{key}' is not present in all sequences.")
            ids = tuple(str(i) for i in range(len(msa)))

    return ids
