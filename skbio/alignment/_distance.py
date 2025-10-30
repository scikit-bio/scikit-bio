# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from typing import Optional, Union, Any, Callable, TYPE_CHECKING
from inspect import isfunction

import numpy as np

from skbio.stats.distance import DistanceMatrix
from skbio.sequence._alphabet import _encode_alphabet

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
    aligned sites and a given distance metric. The resulting distance matrix can be
    used for phylogenetic tree reconstruction (see :mod:`skbio.tree`) or other types of
    distance-based analyses.

    A "site" refers to a position in a sequence that has a valid character from an
    alphabet pre-defined by each metric. Typically, gaps are not considered as sites,
    and will be excluded from calculation. Some metrics are further restricted to
    certain characters, such as the four nucleotides or the 20 basic amino acids, while
    excluding degenerate and/or non-canonical characters.

    Sequences are filtered to aligned sites prior to calculation. This has two working
    modes: By default, alignment positions with at least one invalid character in any
    of the sequences are excluded from all sequences, leaving fully-aligned positions
    across the entire alignment. When ``shared_by_all`` is set to False, this filtering
    is applied to each pair of sequences, and positions with both characters valid are
    retained for calculation, even if there are invalid characters at the same position
    in other sequences.

    Parameters
    ----------
    alignment : TabularMSA of shape (n_sequences, n_positions)
        Multiple sequence alignment.
    metric : str or callable
        The distance metric to apply to each pair of aligned sequences. See
        :mod:`skbio.sequence.distance` for available metrics. Passing metric as a
        string is preferable as this often results in an optimized version of the
        metric being used.
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

    Notes
    -----
    If you wish to apply a custom function between each pair of sequences, use
    ``DistanceMatrix.from_iterable(alignment, function)`` instead.

    This function always strips positions with gap or degenerate characters from the
    aligned sequences. If you intend to compute distances using unstripped sequences,
    consider using ``DistanceMatrix.from_iterable(alignment, metric)``.

    """
    n_seqs = len(alignment)

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
        alphabet = func._alphabet
        func = getattr(skbio.sequence.distance, "_" + metric)

        # func_name = "_" + metric
        # if shared_by_all is False:
        #     func_name += "_pair"
        # func = getattr(skbio.sequence.distance, func_name)

        # Create a 2D matrix of ASCII codes of all sequences.
        # TODO: This can be omitted after optimizing TabularMSA.
        seqs = np.vstack([seq._bytes for seq in alignment])

        # Use all positions (gaps or characters).
        if alphabet is None:
            dm = func(seqs, None, **kwargs)

        # Filter sequences by a given alphabet.
        else:
            valid = _valid_hash(alphabet, alignment.dtype)
            site_mat = valid[seqs]

            # complete deletion
            if shared_by_all:
                site_vec = np.all(site_mat, axis=0)
                seqs = seqs[:, site_vec]
                dm = func(seqs, None, **kwargs)

            # pairwise deletion
            else:
                dm = func(seqs, site_mat, **kwargs)

    # Call a custom function on each pair of sequences (slow).
    elif callable(metric):
        if hasattr(metric, "_alphabet"):
            alphabet = metric._alphabet
        else:
            alphabet = "nogap"

        valid = _valid_hash(alphabet, alignment.dtype)
        site_mat = [valid[seq._bytes] for seq in alignment]

        dm = np.empty((n_seqs * (n_seqs - 1) // 2,))

        # complete deletion
        if shared_by_all:
            site_vec = np.all(site_mat, axis=0)
            seqs = [seq[site_vec] for seq in alignment]
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
                    site_vec = site_mat[i] & site_mat[j]  ## will stripe twice
                    seq1 = alignment[i][site_vec]
                    seq2 = alignment[j][site_vec]
                    dm[pos + j] = metric(seq1, seq2, **kwargs)
                pos += n_seqs - i - 2
    else:
        raise TypeError("`metric` must be a function or a string.")

    return DistanceMatrix(dm, ids=alignment.index.astype(str))


def _valid_hash(alphabet, seqtype):
    """Create a hash table of valid characters for sequence filtering.

    Parameters
    ----------
    alphabet : str, 1D array_like, {'nongap', 'definite', 'canonical'}
        An alphabet of valid characters to be considered by the metric.
    seqtype : class
        Sequence type, which stores grammar information.

    Returns
    -------
    ndarray of bool of shape (128,)
        Boolean mask of ASCII codes of valid characters (True).

    See Also
    --------
    skbio.sequence.distance._metric_specs

    """
    if alphabet in ("nongap", "definite", "canonical"):
        valid = getattr(seqtype, f"_{alphabet}_hash")
    else:
        encoded = _encode_alphabet(alphabet)
        valid = np.zeros((seqtype._num_ascii_codes,), dtype=bool)
        valid[encoded] = True
    return valid


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
