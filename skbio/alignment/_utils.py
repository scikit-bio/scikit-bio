# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.sequence import Sequence, GrammaredSequence, SubstitutionMatrix
from skbio.sequence._alphabet import (
    _encode_alphabet,
    _alphabet_to_hashes,
    _indices_in_observed,
    _indices_in_alphabet,
    _indices_in_alphabet_ascii,
)


# cached identity matrices for alignment with match/mismatch scores
_idmats = {}

# indices of cached identity matrices
_ididxs = {}


def _check_same_type(items):
    """Return the common type of variables.

    Parameters
    ----------
    items : iterable of any
        Variables.

    Returns
    -------
    type
        Common type.

    Raises
    ------
    TypeError
        If there are more than one type.

    """
    # TODO: This function can be moved to skbio.util._misc.
    dtype = type(items[0])
    for item in items[1:]:
        if type(item) is not dtype:
            raise TypeError("Variables are of different types.")
    return dtype


def _parse_seqs(seqs):
    """Parse sequences for alignment.

    Parameters
    ----------
    seqs : iterable of iterable of scalar
        Sequences.

    Returns
    -------
    list of ndarray of uint8 or int
        Encoded sequences.
    type, int, or "ascii"
        Sequence type | alphabet is ASCII | alphabet size

    Notes
    -----
    It converts sequences into arrays of indices in a hypothetical subsitution matrix,
    and suggests how to construct this subsitution matrix. Three scenarios:

        1) GrammaredSequence, with a finite alphabet.
           => size of alphabet
        2) Non-grammared Sequence, characters and integers within ASCII range (0-127).
           => size = 128
        3) Unicode characters, arbitrary numbers, hashables.
           => size = count of unique elements

    """
    seqtype = _check_same_type(seqs)
    if issubclass(seqtype, GrammaredSequence):
        return [x._bytes for x in seqs], seqtype
    elif issubclass(seqtype, Sequence):
        return [x._bytes for x in seqs], "ascii"
    try:
        seqs = [_encode_alphabet(x) for x in seqs]
    except (TypeError, ValueError, UnicodeEncodeError):
        seqs, uniq = _indices_in_observed(seqs)
        return seqs, uniq.size
    else:
        return seqs, "ascii"


def _parse_seqs_submat(seqs, submat):
    """Parse sequences with substitution matrix.

    Parameters
    ----------
    seqs : iterable of iterable of scalar
        Sequences.
    submat : SubstitutionMatrix
        Substitution matrix.

    Returns
    -------
    list of ndarray of uint8 or int
        Indices of sequence characters in the substitution matrix.

    """
    seqtype = _check_same_type(seqs)
    if issubclass(seqtype, Sequence):
        return [x.to_indices(submat, mask_gaps=False) for x in seqs]
    if submat._is_ascii:
        try:
            seqs = [_encode_alphabet(x) for x in seqs]
        except (TypeError, ValueError, UnicodeEncodeError):
            raise ValueError("Cannot encode ASCII characters.")
        return [_indices_in_alphabet_ascii(x, submat._char_hash) for x in seqs]
    else:
        return [_indices_in_alphabet(x, submat._char_map) for x in seqs]


def _prep_seqs_submat(seqs, sub_score, dtype=np.float64):
    """Prepare sequences and substitution matrix for alignment.

    Parameters
    ----------
    seqs : iterable of Sequence or sequence
        Sequences.
    sub_score : tuple of (float, float), SubstitutionMatrix, or str
        Score of a substitution. May be two numbers (match, mismatch), a substitution
        matrix, or its name.
    dtype : type, optional
        Default floating-point data type (np.float32 or np.float64) to use if it cannot
        be automatically inferred.

    Returns
    -------
    list of ndarray of uint8 or int
        Indices of sequence characters in the substitution matrix.
    ndarray of float of shape (n_alphabet, n_alphabet)
        Substitution matrix.

    """
    if isinstance(sub_score, str):
        sub_score = SubstitutionMatrix.by_name(sub_score)
    if isinstance(sub_score, SubstitutionMatrix):
        seqs = _parse_seqs_submat(seqs, sub_score)
        submat = sub_score._data
    else:
        match, mismatch = sub_score
        seqs, key = _parse_seqs(seqs)

        # create a temporary matrix of arbitrary size
        if isinstance(key, int):
            submat = np.full((key, key), mismatch, dtype=dtype)
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
                alphabet = [chr(i) for i in range(128)]
                seqs = [x.astype(np.intp) for x in seqs]
            else:
                alphabet = sorted(key.alphabet)
                # create index and translate
                idx = _ididxs[key] = _alphabet_to_hashes(alphabet)
                seqs = [idx[x] for x in seqs]

            n = len(alphabet)
            submat = np.full((n, n), mismatch, dtype=dtype)
            np.fill_diagonal(submat, match)
            _idmats[key] = submat

    return seqs, submat


def _prep_gap_cost(gap_cost, dtype=np.float64):
    """Prepare gap penalty for alignment.

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

    """
    if np.isscalar(gap_cost):
        gap_open, gap_extend = 0, gap_cost
    else:
        gap_open, gap_extend = gap_cost
    return dtype(gap_open), dtype(gap_extend)
