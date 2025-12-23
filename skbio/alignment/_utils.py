# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from collections.abc import Iterable
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

from skbio.alignment import TabularMSA, AlignPath
from skbio.sequence import Sequence, GrammaredSequence, SubstitutionMatrix
from skbio.sequence._alphabet import (
    _encode_alphabet,
    _alphabet_to_hashes,
    _indices_in_observed,
)


# This could be exposed as a public API.
SequenceLike: TypeAlias = (
    Sequence | str | bytes | Iterable[str | bytes | int | float | bool] | NDArray
)


# This could be exposed as a public API.
AlignmentLike: TypeAlias = (
    TabularMSA | Iterable[SequenceLike] | tuple[AlignPath, Iterable[SequenceLike]]
)


# cached identity matrices for alignment with match/mismatch scores
_idmats: dict[str | type, NDArray] = {}

# indices of cached identity matrices
_ididxs: dict[type, NDArray] = {}


def encode_sequences(
    seqs: Iterable[SequenceLike],
    sub_score: tuple[float, float] | SubstitutionMatrix | str,
    aligned: bool = False,
    dtype: type = np.float32,
    gap_chars: str = "-",
    not_empty: bool = True,
) -> tuple[list[NDArray], NDArray, NDArray]:
    """Encode sequences for alignment operations.

    This function transforms sequences into indices in a substitution matrix to
    facilitate subsequent alignment operations.

    Parameters
    ----------
    seqs : iterable of Sequence, str, or sequence of scalar
        Input sequences.
    sub_score : tuple of (float, float), SubstitutionMatrix, or str
        Substitution scoring method. Can be two numbers (match, mismatch), a
        substitution matrix, or its name.
    aligned : bool, optional
        Whether sequences are aligned. If True, sequences will be checked for equal
        length and gaps will be encoded into a Boolean mask.
    dtype : type, optional
        Default floating-point data type (np.float32 or np.float64) to use if it cannot
        be automatically inferred.
    gap_chars : str, optional
        Characters that should be treated as gaps in aligned sequences. Default is "-".
    not_empty : bool, optional
        If True (default), each sequence must not be empty.

    Returns
    -------
    seqs : list of ndarray of intp
        Indices of sequence characters in the substitution matrix.
    submat : ndarray of float of shape (n_alphabet, n_alphabet)
        Substitution matrix.
    gaps : ndarray of bool of shape (n_sequences, n_positions), optional
        Boolean mask of gap positions.

    Raises
    ------
    TypeError
        If sequences are of heterogeneous types.
    ValueError
        If there is no sequence.
    ValueError
        If any sequence has a length of zero.

    See Also
    --------
    encode_alignment

    Notes
    -----
    This function is currently private but it could be exposed as a public API.

    """
    # Determine type of sequences. They can be skbio sequences (grammared or not),
    # raw strings or any iterables of scalars. Heterogeneous sequences are not allowed.
    seqtype = _check_seqtype(seqs)
    grammared = issubclass(seqtype, GrammaredSequence)

    # Determine substitution scoring method. It can be a specified substitution matrix
    # or match & mismatch scores. For the later scenario, an identity matrix will be
    # later constructed or retrieved from the cache.
    if isinstance(sub_score, str):
        submat = SubstitutionMatrix.by_name(sub_score)
        is_submat = True
    elif isinstance(sub_score, SubstitutionMatrix):
        submat = sub_score
        is_submat = True
    else:
        is_submat = False
        match, mismatch = sub_score
        match, mismatch = dtype(match), dtype(mismatch)

    # Determine whether sequences should be encoded as ASCII codes, which skbio has
    # optimizations for.
    # This requires that 1) sequences are ASCII, and 2) substitution matrix (if
    # provided) supports ASCII.
    # If yes, convert sequences into ASCII codes. Otherwise, keep sequences as strings
    # or whatever iterables in the original form.
    if not is_submat or submat.is_ascii:
        if issubclass(seqtype, Sequence):
            is_ascii = True
            seqs = [x._bytes for x in seqs]
        else:
            try:
                seqs = [_encode_alphabet(x) for x in seqs]
            except (TypeError, ValueError, UnicodeEncodeError):
                if is_submat:
                    raise ValueError(
                        "Substitution matrix has an ASCII alphabet, but sequences "
                        "cannot be fully encoded into ASCII codes."
                    )
                is_ascii = False
            else:
                is_ascii = True
    else:
        is_ascii = False
        if issubclass(seqtype, Sequence):
            seqs = [str(x) for x in seqs]

    # Identify gaps in aligned sequences. The output is a 2D Boolean mask representing
    # gap positions (0: character, 1: gap). This process also validates that all
    # sequences in the alignment are of equal length.
    gaps = None
    if aligned:
        if is_ascii:
            if grammared:
                gap_codes = seqtype._gap_codes
            else:
                gap_codes = [ord(x) for x in gap_chars]
        else:
            if grammared:
                gap_codes = list(seqtype.gap_chars)
            else:
                gap_codes = list(gap_chars)

        gaps = _mask_gaps(seqs, gap_codes)

    # Index sequences according to a given substitution matrix.
    if is_submat:
        # Determine wildcard if available.
        wild = None
        if grammared:
            wild = getattr(seqtype, "wildcard_char", None)
            if is_ascii and wild is not None:
                wild = ord(wild)

        if is_ascii:
            seqs = _map_chars_ascii(seqs, submat._char_hash, wild=wild)
        else:
            seqs = _map_chars(seqs, submat._char_map, wild=wild)

        # Validate that all non-gap positions have valid indices.
        _check_indices(seqs, gaps)

        submat = submat._data

    # Construct a substitution matrix based on match/mismatch scores, or retrieve one
    # from cache.
    else:
        if is_ascii:
            key = seqtype if grammared else "ascii"
        else:
            seqs, uniq = _indices_in_observed(seqs)
            key = uniq.size
        seqs, submat = prep_identity_matrix(seqs, key, match, mismatch)

    if not_empty:
        for i, seq in enumerate(seqs):
            if seq.size == 0:
                raise ValueError(f"Sequence {i + 1} has a length of zero.")

    return seqs, submat, gaps


def encode_alignment(
    aln: AlignmentLike,
    sub_score: tuple[float, float] | SubstitutionMatrix | str,
    gap_chars: str = "-",
) -> tuple[list[NDArray], NDArray, NDArray, NDArray]:
    """Encode sequences for alignment operations.

    This function transforms an alignment into a 2D array of indices in a
    substitution matrix, and run-length encoded gaps, which will facilitate
    subsequent alignment operations.

    Parameters
    ----------
    aln : TabularMSA, iterable, or (AlignPath, iterable)
        Input alignment.
    sub_score : tuple of (float, float), SubstitutionMatrix, or str
        Substitution scoring method. Can be two numbers (match, mismatch), a
        substitution matrix, or its name.
    gap_chars : str, optional
        Characters that should be treated as gaps in aligned sequences.

    Returns
    -------
    seqs : list of ndarray of intp
        Indices of sequence characters in the substitution matrix.
    submat : ndarray of float of shape (n_alphabet, n_alphabet)
        Substitution matrix.
    bits : ndarray of bool of shape (n_sequences, n_segments)
        Gap status of segments in the alignment path.
    lens : ndarray of int of shape (n_segments,)
        Lengths of segments in the alignment path.

    Raises
    ------
    TypeError
        If sequences are of heterogeneous types.
    ValueError
        If the alignment has no sequence.
    ValueError
        If the alignment has no column (position).

    See Also
    --------
    encode_sequences

    Notes
    -----
    This function is currently private but it could be exposed as a public API.

    """
    # There are two scenarios:
    is_path = False
    if isinstance(aln, (list, tuple)) and len(aln) == 2:
        path, seqs = aln
        if isinstance(path, AlignPath) and isinstance(seqs, Iterable):
            is_path = True

    # 1. Input is an alignment path and original (unaligned) sequences.
    if is_path:
        seqs, submat, _ = encode_sequences(
            seqs, sub_score, aligned=False, gap_chars=gap_chars
        )
        seqs, _, bits, lens = path._to_matrices(seqs)

    # 2. Input is aligned sequences (with gaps inside them).
    else:
        seqs, submat, gaps = encode_sequences(
            aln, sub_score, aligned=True, gap_chars=gap_chars, not_empty=False
        )
        seqs = np.vstack(seqs)
        if seqs.shape[1] == 0:
            raise ValueError("The alignment has a length of zero.")
        bits, lens = _get_align_path(gaps)

    return seqs, submat, bits, lens


def prep_gapcost(gap_cost, dtype=np.float32):
    """Prepare gap penalty method for alignment.

    Parameters
    ----------
    gap_cost : float or tuple of (float, float), optional
        Penalty of a gap. May be one (linear) or two numbers (affine). See
        :func:`align_score` for instructions and rationales. Default is -2.
    dtype : type, optional
        Floating-point data type (np.float32 or np.float64) to use.

    Returns
    -------
    float
        Gap opening penalty.
    float
        Gap extension penalty.

    Notes
    -----
    Gap penalties are usually positive, representing subtractions from the alignment
    score. Although scikit-bio does not prohibit negative gap penalties, it should be
    noted that negative penalties (positive scores) might lead to all-gap alignments,
    especially when gaps are scored higher than substitutions.

    """
    if np.isscalar(gap_cost):
        gap_open, gap_extend = 0, gap_cost
    else:
        gap_open, gap_extend = gap_cost
    return dtype(gap_open), dtype(gap_extend)


def prep_identity_matrix(seqs, key, match, mismatch):
    """Prepare an identity matrix based on match/mismatch scores.

    Parameters
    ----------
    seqs : list of ndarray of uint8 or int
        Encoded sequences.
    key: type, int, or "ascii"
        Sequence type | alphabet is ASCII | alphabet size
    match : float
        Match score.
    mismatch : float
        Mismatch score.

    Returns
    -------
    list of ndarray of intp
        Indices of sequence characters in the substitution matrix.
    ndarray of float of shape (n_alphabet, n_alphabet)
        Substitution matrix.

    Notes
    -----
    An identity matrix is created on the first use and saved in the cache for reuse.
    There are three scenarios:

    1. Each GrammaredSequence type (such as DNA and Protein) will have its own identity
       matrix that covers its specific alphabet.
    2. For non-grammared Sequence, string and other types that can be ASCII-encoded, a
       single identity matrix with all 128 ASCII codes will be created and cached.
    3. For other types, an identity matrix will be created each time and not cached.
       It will cover all unique elements in the sequences.

    """
    # create a temporary matrix of arbitrary size
    if isinstance(key, int):
        submat = np.full((key, key), mismatch)
        np.fill_diagonal(submat, match)

    # check existing matrix in cache and update it if needed
    elif key in _idmats:
        submat = _idmats[key]
        if submat[0, 1] != mismatch:
            submat[:] = mismatch
        if submat[0, 0] != match:
            np.fill_diagonal(submat, match)
        if key == "ascii":
            # cast ASCII codes (uint8) into indices (intp)
            seqs = [x.astype(np.intp) for x in seqs]
        else:
            # translate ASCII codes to indices in matrix
            idx = _ididxs[key]
            seqs = [idx[x] for x in seqs]

    # create a new matrix and save it to cache
    else:
        if key == "ascii":
            # TODO: This can be further optimized in `SubstitutionMatrix`.
            alphabet = [chr(i) for i in range(128)]
            seqs = [x.astype(np.intp) for x in seqs]
        else:
            alphabet = sorted(key.alphabet)
            # create index and translate
            idx = _ididxs[key] = _alphabet_to_hashes(alphabet)
            seqs = [idx[x] for x in seqs]

        n = len(alphabet)
        submat = np.full((n, n), mismatch)
        np.fill_diagonal(submat, match)
        _idmats[key] = submat

    return seqs, submat


def _check_seqtype(seqs):
    """Check if sequences are homogeneous, and return the common type."""
    if isinstance(seqs, TabularMSA):
        dtype = seqs.dtype
        if seqs.shape[0] == 0:
            raise ValueError("No sequence is provided.")
        return dtype
    try:
        dtype = type(seqs[0])
    except IndexError:
        raise ValueError("No sequence is provided.")
    for seq in seqs[1:]:
        if type(seq) is not dtype:
            raise TypeError("Sequences are of different types.")
    return dtype


def _mask_gaps(seqs, gap_codes):
    """Mask gap positions in aligned sequences."""
    if isinstance(seqs[0], str):
        gaps = [np.isin(list(x), gap_codes) for x in seqs]
    else:
        gaps = [np.isin(x, gap_codes) for x in seqs]
    try:
        return np.vstack(gaps)
    except ValueError:
        raise ValueError("Sequence lengths do not match.")


def _map_chars_ascii(seqs, mapping, wild=None):
    """Map characters in sequences to indices.

    `mapping` is an array of ASCII codes (n=128).

    This function can be considered as a refined and batched version of
    `skbio.sequence._alphabet._indices_in_alphabet_ascii`.

    """
    seqs = [mapping[x] for x in seqs]
    if wild is not None and (wild := mapping[wild]) != -1:
        for seq in seqs:
            seq[seq == -1] = wild
    return seqs


def _map_chars(seqs, mapping, wild=None):
    """Map characters in sequences to indices.

    `mapping` is a dictionary of original characters.

    This function can be considered as a refined and batched version of
    `skbio.sequence._alphabet._indices_in_alphabet`.

    """
    wild = -1 if wild is None else mapping.get(wild, -1)
    return [np.array([mapping.get(x, wild) for x in y], dtype=np.intp) for y in seqs]


def _check_indices(seqs, gaps=None):
    """Validate that all non-gap positions have valid indices.

    `gaps` is a 2D Boolean mask whereas `seqs` are a list of equal-size arrays (they
    will be vstack'ed later).

    """
    msg = "Sequence {} contain character(s) absent from the substitution matrix."
    if gaps is None:
        for i, seq in enumerate(seqs):
            if (seq == -1).any():
                raise ValueError(msg.format(i + 1))
    else:
        for i, (seq, gap) in enumerate(zip(seqs, gaps)):
            if (seq[~gap] == -1).any():
                raise ValueError(msg.format(i + 1))


def _get_align_path(bits):
    """Calculate the path of an alignment.

    This function is similar to `AlignPath.from_tabular`, except that bits don't need
    to be packed into uint8's.

    """
    idx = np.append(0, np.where((bits[:, :-1] != bits[:, 1:]).any(axis=0))[0] + 1)
    lens = np.append(idx[1:] - idx[:-1], bits.shape[1] - idx[-1])
    bits = bits[:, idx]
    return bits, lens
