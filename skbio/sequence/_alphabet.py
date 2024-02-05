# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np


def _encode_alphabet(alphabet):
    """Encode an alphabet as a vector of ASCII code points.

    Parameters
    ----------
    alphabet : str, list, tuple or 1D np.ndarray
        Input alphabet. Must consist of single ASCII characters. Elements may
        be string or byte characters, or integers representing code points.

    Returns
    -------
    1D np.ndarray of np.uint8
        Vector of ASCII code points representing the alphabet.

    Raises
    ------
    TypeError
        If alphabet or its components are of a wrong data type.
    ValueError
        If some elements are not single characters.
    ValueError
        If some code points are beyond the ASCII range.
    UnicodeEncodeError
        If some characters are beyond the ASCII range.

    Notes
    -----
    ASCII has 128 code points (0 to 127) [1]_ (not to be confused with extended
    ASCII). Therefore, output values are within the range of [0, 127].

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/ASCII

    """
    errmsg = "Alphabet is of an invalid data type."

    # string
    if isinstance(alphabet, str):
        alphabet = alphabet.encode("ascii")
        return np.frombuffer(alphabet, dtype=np.uint8)

    # list or tuple
    elif isinstance(alphabet, (list, tuple)):
        alphabet = np.array(alphabet)

    # 1d numpy array
    elif not isinstance(alphabet, np.ndarray):
        raise TypeError(errmsg)
    if alphabet.ndim != 1:
        raise TypeError(errmsg)
    dtype = alphabet.dtype

    # integers represent ascii code points
    if np.issubdtype(dtype, np.integer):
        # ascii code points are within [0, 127]
        if np.all((alphabet >= 0) & (alphabet <= 127)):
            if dtype is np.uint8:
                return alphabet

            # cast data type to uint8
            else:
                return alphabet.astype(np.uint8)
        else:
            raise ValueError("Not all code points are within the ASCII range.")

    # encode strings as ascii characters
    elif np.issubdtype(dtype, np.str_):
        alphabet = np.char.encode(alphabet, encoding="ascii")

    # bytes are already encoded
    elif not np.issubdtype(dtype, np.bytes_):
        raise TypeError(errmsg)

    # must be single characters
    if not (np.char.str_len(alphabet) == 1).all():
        raise ValueError("Not all elements are single characters.")
    return alphabet.view(np.uint8)


def _alphabet_to_hashes(alphabet):
    """Convert an alphabet into a hash table of ASCII code points to indices.

    Parameters
    ----------
    alphabet : iterable
        Input alphabet. Must consist of single ASCII characters.

    Returns
    -------
    np.ndarray of np.uint8 of shape (128,)
        Hash table of ASCII code points to indices.

    Raises
    ------
    ValueError
        If the absence character is not in the alphabet.
    ValueError
        If one or multiple characters in the sequence are absent from the
        alphabet, whereas `absence` is not set.

    See Also
    --------
    _indices_in_alphabet_ascii

    Notes
    -----
    The resulting data structure enables efficient conversion of a sequence
    into indices of characters in an alphabet.

    The hash table has a constant size of 128, which is the total number of
    ASCII characters.

    Code points absent from the alphabet are filled with 255, which is beyond
    the range of ASCII characters, hence the maximum index in the alphabet.

    """
    idx = _encode_alphabet(alphabet)
    res = np.full(128, 255, dtype=np.uint8)
    res[idx] = np.arange(idx.size)
    return res


def _indices_in_alphabet(seq, alphabet, wildcard=None):
    """Convert a sequence into indices of characters in an alphabet.

    Parameters
    ----------
    seq : iterable
        Input sequence.
    alphabet : dict or iterable
        Input alphabet. Can be a dictionary of characters to indices, or an
        iterable of other types from which the dictionary will be constructed.
    wildcard : hashable, optional
        Character to replace any characters that are absent from the alphabet.
        If omitted, will raise an error if the latter characters exist.

    Returns
    -------
    1D np.ndarray of int
        Vector of indices of characters in an alphabet.

    Raises
    ------
    ValueError
        If the wildcard character is not in the alphabet.
    ValueError
        If one or multiple characters in the sequence are absent from the
        alphabet, whereas `wildcard` is not set.

    See Also
    --------
    _indices_in_alphabet_ascii

    Notes
    -----
    This function is versatile to the type of characters.

    """
    if not isinstance(alphabet, dict):
        alphabet = {x: i for i, x in enumerate(alphabet)}
    pos = list(map(alphabet.get, seq))
    if wildcard is not None:
        try:
            wildcard = alphabet[wildcard]
        except KeyError:
            raise ValueError(
                f'Wildcard character "{wildcard}" is not in the ' "alphabet."
            )
        pos = [wildcard if x is None else x for x in pos]
    elif None in pos:
        raise ValueError(
            "One or multiple characters in the sequence are "
            "absent from the alphabet."
        )
    return np.array(pos)


def _indices_in_alphabet_ascii(seq, alphabet, wildcard=None):
    """Convert a sequence into indices of characters in an ASCII alphabet.

    Parameters
    ----------
    seq : 1D np.ndarray of int
        Input sequence as ASCII code points.
    alphabet : np.ndarray of shape (128,) of int
        Input alphabet as a hash table of all ASCII code points to character
        indices, or 255 if absent from the alphabet.
    wildcard : int, optional
        Code point of character to replace any characters that are absent from
        the alphabet. If omitted, will raise an error if such characters exist.

    Returns
    -------
    1D np.ndarray of uint8
        Vector of indices of characters in an alphabet.

    Raises
    ------
    ValueError
        If the wildcard character is not in the alphabet.
    ValueError
        If one or multiple characters in the sequence are absent from the
        alphabet, whereas `wildcard` is not set.

    See Also
    --------
    _indices_in_alphabet
    _alphabet_to_hashes

    Notes
    -----
    This function is optimized for single ASCII characters.

    """
    pos = alphabet[seq]
    absent = pos == 255
    if absent.any():
        if wildcard is None:
            raise ValueError(
                "One or multiple characters in the sequence are "
                "absent from the alphabet."
            )
        try:
            assert (wild := alphabet[wildcard]) != 255
        except AssertionError:
            raise ValueError(
                f'Wildcard character "{chr(wildcard)}" is not in ' "the alphabet."
            )
        pos = np.where(absent, wild, pos)
    return pos


def _indices_in_observed(seqs):
    """Convert sequences into vectors of indices in observed characters.

    Parameters
    ----------
    seqs : iterable of iterable
        Input sequences.

    Returns
    -------
    list of 1D np.ndarray
        Vectors of indices representing the sequences.
    1D np.ndarray
        Sorted vector of unique characters observed in the sequences.

    """
    # This function uses np.unique to extract unique characters and their
    # indices. It applies np.unique on individual sequences, then merges
    # results. This design is to avoid concatenating too many sequences.
    alpha_lst, index_lst = zip(
        *[
            np.unique(tuple(x) if isinstance(x, str) else x, return_inverse=True)
            for x in seqs
        ]
    )
    alpha_union, index_union = np.unique(np.concatenate(alpha_lst), return_inverse=True)
    index_bounds = np.cumsum([x.size for x in alpha_lst])[:-1]
    index_chunks = np.split(index_union, index_bounds)
    index_lst_trans = [x[y] for x, y in zip(index_chunks, index_lst)]
    return index_lst_trans, alpha_union
