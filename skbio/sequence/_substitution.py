# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.util._decorator import classonlymethod
from skbio.stats.distance import DissimilarityMatrix
from skbio.sequence._alphabet import _alphabet_to_hashes


class SubstitutionMatrix(DissimilarityMatrix):
    """Scoring matrix between characters in biological sequences.

    Parameters
    ----------
    alphabet : iterable
        Characters that constitute the alphabet.
    scores : 2D array-like
        Scores of substitutions from one character (row, or axis=0) to another
        character (column, or axis=1).
    kwargs : dict
        Additional arguments for the ``DissimilarityMatrix`` constructor.

    See Also
    --------
    skbio.stats.distance.DissimilarityMatrix

    Notes
    -----
    A substitution matrix (a.k.a. replacement matrix) scores the substitution
    of each character by each other character and itself in an alphabet. The
    score usually represents the rate of substitution over evolutionary time in
    a biological sequence. A higher score usually indicates higher similarity
    in chemical properties or functional roles of two molecules, therefore a
    mutation from one to the other is easier. In sequence alignment, the score
    can measure the likelihood that a pair of aligned characters are homologous
    rather than by chance.

    This class provides a generalized interface for substitution matrices. The
    alphabet usually consists of individual characters, such as nucleotides or
    amino acids, but it can be generalized to any iterable of scalars (numbers,
    strings, etc.). Therefore, you may use this class to construct substitution
    matrices of complicated biological units (such as codons or non-canonical
    amino acids). The score matrix may be symmetric, as many existing matrices
    are, or asymmetric, where the score of one character substituted by another
    is unequal to the other way around. Only square matrices (i.e., numbers of
    rows and columns are equal) are supported.

    Multiple commonly used nucleotide and amino acid substitution matrices are
    pre-defined and can be referred to by name. Examples include NUC.4.4 for
    nucleotides, and variants of BLOSUM and PAM matrices for amino acids.

    ``SubstitutionMatrix`` is a subclass of ``DissimilarityMatrix``. Therefore,
    all attributes and methods of the latter also apply to the former.

    Examples
    --------
    >>> from skbio import SubstitutionMatrix
    >>> mat = SubstitutionMatrix('ACGT', np.array([
    ...     [2, -1, -1, -1],
    ...     [-1, 2, -1, -1],
    ...     [-1, -1, 2, -1],
    ...     [-1, -1, -1, 2]]))
    >>> mat.alphabet
    ('A', 'C', 'G', 'T')
    >>> mat.scores
    array([[ 2., -1., -1., -1.],
           [-1.,  2., -1., -1.],
           [-1., -1.,  2., -1.],
           [-1., -1., -1.,  2.]])
    >>> mat['A', 'T']
    -1.0
    >>> mat['G', 'G']
    2.0
    >>> blosum62 = SubstitutionMatrix.by_name('BLOSUM62')

    """

    @property
    def alphabet(self):
        """Alphabet of the substitution matrix.

        Each element (character) corresponds to one row/column in the matrix.

        Returns
        -------
        tuple
            Alphabet of the substitution matrix.

        Notes
        -----
        This is an alias of ``ids``.

        """
        return self._ids

    @property
    def scores(self):
        """Matrix of substitution scores.

        Each value corresponds to the score of substituting the row character
        with the column character.

        Returns
        -------
        2D np.ndarray
            Matrix of substitution scores.

        Notes
        -----
        This is an alias of ``data``.

        """
        return self._data

    @property
    def is_ascii(self):
        """Whether alphabet consists of single ASCII characters.

        `True` if every character in the alphabet can be represented by a
        single ASCII character within code point range 0 to 255.

        Returns
        -------
        bool
            Whether alphabet consists of single ASCII characters.

        """
        return self._is_ascii

    def __init__(self, alphabet, scores, **kwargs):
        """Initialize a substitution matrix object."""
        super().__init__(scores, alphabet, **kwargs)

        # `_char_map`: dictionary of characters to indices in the alphabet.
        # It is to enable efficient conversion of sequences into indices.
        # It is compatible with generalized sequences.
        self._char_map = {x: i for i, x in enumerate(alphabet)}

        # `_is_ascii`: whether alphabet can be encoded as ASCII characters.
        # `_char_hash`: hash table of ASCII code points to indices in the
        # alphabet. It is to enable efficient conversion of sequences into
        # into indices. It is optimized for ASCII sequences.
        try:
            hash_ = _alphabet_to_hashes(alphabet)
        except (TypeError, ValueError, UnicodeEncodeError):
            self._is_ascii = False
        else:
            self._is_ascii = True
            self._char_hash = hash_

    def to_dict(self):
        """Create a 2D dictionary from the substitution matrix.

        Returns
        -------
        dict of dict
            2D dictionary constructed from the substitution matrix.

        """
        return {
            id_: dict(zip(self._ids, row)) for id_, row in zip(self._ids, self._data)
        }

    @classonlymethod
    def from_dict(cls, dictionary):
        """Create a substitution matrix from a 2D dictionary.

        Parameters
        ----------
        dictionary : dict of dict
            2D dictionary of substitution scores from outer characters to inner
            characters.

        Returns
        -------
        SubstitutionMatrix
            Substitution matrix constructed from the dictionary.

        Raises
        ------
        ValueError
            If outer and inner characters are inconsistent.
        ValueError
            If scores are not numbers.

        Examples
        --------
        >>> from skbio import SubstitutionMatrix
        >>> d = {'a': {'a': 1, 'b': 0, 'c': 0},
        ...      'b': {'a': 0, 'b': 1, 'c': 0},
        ...      'c': {'a': 0, 'b': 0, 'c': 1}}
        >>> mat = SubstitutionMatrix.from_dict(d)
        >>> mat.alphabet
        ('a', 'b', 'c')
        >>> mat.scores
        array([[ 1.,  0.,  0.],
               [ 0.,  1.,  0.],
               [ 0.,  0.,  1.]])

        """
        alphabet, rows = zip(*dictionary.items())
        alphabet_set = set(alphabet)
        idmap = {x: i for i, x in enumerate(alphabet)}
        scores = np.zeros((n := len(alphabet), n))
        for i, row in enumerate(rows):
            if set(row) != alphabet_set:
                raise ValueError(
                    "The outer and inner layers of the dictionary"
                    " must have the same set of keys."
                )
            for key, value in row.items():
                scores[i][idmap[key]] = float(value)

        return cls(alphabet, scores)

    @classonlymethod
    def identity(cls, alphabet, match, mismatch):
        """Create an identity substitution matrix.

        All matches and mismatches will have the identical scores,
        respectively, regardless of the character.

        Parameters
        ----------
        alphabet : iterable
            Characters that constitute the alphabet.
        match : int or float
            Score assigned to all matches.
        mismatch : int or float
            Score assigned to all mismatches.

        Returns
        -------
        SubstitutionMatrix
            Substitution matrix constructed from the alphabet and scores.

        Examples
        --------
        >>> from skbio import SubstitutionMatrix
        >>> mat = SubstitutionMatrix.identity('ACGT', 1, -2)
        >>> mat.alphabet
        ('A', 'C', 'G', 'T')
        >>> mat.scores
        array([[ 1., -2., -2., -2.],
               [-2.,  1., -2., -2.],
               [-2., -2.,  1., -2.],
               [-2., -2., -2.,  1.]])

        """
        alphabet = tuple(alphabet)
        scores = np.identity(len(alphabet)) * (match - mismatch) + mismatch
        return cls(alphabet, scores)

    @classonlymethod
    def by_name(cls, name):
        """Load a pre-defined substitution matrix by its name.

        Parameters
        ----------
        name : str
            Name of the substitution matrix.

        Returns
        -------
        SubstitutionMatrix
            Named substitution matrix.

        Raises
        ------
        ValueError
            If named substitution matrix does not exist.

        See Also
        --------
        get_names

        Notes
        -----
        Names are case-insensitive. For instance, `BLOSUM62` and `blosum62`
        point to the same substitution matrix.

        Available substitution matrix names can be obtained by ``get_names``.
        Currently, the following names are supported:

        - `NUC.4.4` (a.k.a. DNAfull): A nucleotide substitution matrix covering
          all definite and degenerate nucleotides.

        - Point Accepted Mutation (PAM) [1]_: A set of amino acid substitution
          matrices, including `PAM30`, `PAM70` and `PAM250`.

        - BLOcks SUbstitution Matrix (BLOSUM) [2]_: A set of amino acid
          substitution matrices, including `BLOSUM45`, `BLOSUM50`, `BLOSUM62`,
          `BLOSUM80` and `BLOSUM90`.

        References
        ----------
        .. [1] Dayhoff, M., Schwartz, R., & Orcutt, B. (1978). A model of
           evolutionary change in proteins. Atlas of protein sequence and
           structure, 5, 345-352.
        .. [2] Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution
           matrices from protein blocks. Proceedings of the National Academy of
           Sciences, 89(22), 10915-10919.

        Examples
        --------
        >>> from skbio import SubstitutionMatrix
        >>> mat = SubstitutionMatrix.by_name('BLOSUM62')
        >>> len(mat.alphabet)
        24
        >>> mat['M', 'K']
        -1.0

        """
        try:
            return _named_substitution_matrices[name]
        except KeyError:
            name_lower = name.lower()
            for key, value in _named_substitution_matrices.items():
                if name_lower == key.lower():
                    return value
            raise ValueError(f'Substitution matrix "{name}" does not exist.')

    @classonlymethod
    def get_names(cls):
        """List names of pre-defined substitution matrices.

        Returns
        -------
        list of str
            Names of pre-defined substitution matrices.

        See Also
        --------
        by_name

        """
        return list(_named_substitution_matrices.keys())


def _matrix_to_vector(mat):
    """Flatten a square matrix to a vector of the upper triangle and diagonal."""
    assert len(mat.shape) == 2
    assert mat.shape[0] == mat.shape[1]
    return mat[np.triu_indices(len(mat))]


def _vector_to_matrix(vec):
    """Revert a vector representing a flattened matrix to square form."""
    assert len(vec.shape) == 1
    # a square matrix of shape (n, n) will have n * (n + 1) / 2 elements in the
    # flattened vector; the following code reverses this equation to obtain the
    # original shape of the matrix
    n = (np.sqrt(1 + 8 * len(vec)) - 1) / 2
    assert n == (n := int(n))
    mat = np.zeros((n, n))
    mat[np.triu_indices(n)] = vec
    return mat + np.triu(mat, k=1).T


# Defined according to the matrices hosted at the NCBI FTP server:
# https://ftp.ncbi.nlm.nih.gov/blast/matrices/
# fmt: off
_named_substitution_matrices = {
    # NUC.4.4, a.k.a. DNAfull
    "NUC.4.4": SubstitutionMatrix(
        "ATGCSWRYKMBVHDN",
        _vector_to_matrix(
            np.array(
                [5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2, 5, -4, -4,
                -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, 5, -4, 1, -4, 1, -4, 1,
                -4, -1, -1, -4, -1, -2, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2,
                -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1, -1, -2, -2, -2, -2, -3,
                -3, -1, -1, -1, -1, -4, -2, -2, -3, -1, -3, -1, -1, -1, -2, -2, -1,
                -3, -1, -3, -1, -1, -4, -1, -3, -3, -1, -1, -1, -3, -1, -1, -3, -1,
                -1, -2, -2, -2, -1, -1, -2, -2, -1, -1, -2, -1, -1, -1, -1]
            )
        ),
        validate=False,
    ),
    # Point Accepted Mutation (PAM)
    "PAM30": SubstitutionMatrix(
        "ARNDCQEGHILKMFPSTWYVBZX*",
        _vector_to_matrix(
            np.array(
                [ 6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5, -8, -2, 0, -1,
                -13, -8, -2, -3, -3, -3, -17, 8, -6, -10, -8, -2, -9, -9, -2, -5,
                -8, 0, -4, -9, -4, -3, -6, -2, -10, -8, -7, -4, -6, -17, 8, 2, -11,
                -3, -2, -3, 0, -5, -7, -1, -9, -9, -6, 0, -2, -8, -4, -8, 6, -3,
                -3, -17, 8, -14, -2, 2, -3, -4, -7, -12, -4, -11, -15, -8, -4, -5,
                -15, -11, -8, 6, 1, -5, -17, 10, -14, -14, -9, -7, -6, -15, -14,
                -13, -13, -8, -3, -8, -15, -4, -6, -12, -14, -9, -17, 8, 1, -7, 1,
                -8, -5, -3, -4, -13, -3, -5, -5, -13, -12, -7, -3, 6, -5, -17, 8,
                -4, -5, -5, -9, -4, -7, -14, -5, -4, -6, -17, -8, -6, 1, 6, -5,
                -17, 6, -9, -11, -10, -7, -8, -9, -6, -2, -6, -15, -14, -5, -3, -5,
                -5, -17, 9, -9, -6, -6, -10, -6, -4, -6, -7, -7, -3, -6, -1, -1,
                -5, -17, 8, -1, -6, -1, -2, -8, -7, -2, -14, -6, 2, -6, -6, -5,
                -17, 7, -8, 1, -3, -7, -8, -7, -6, -7, -2, -9, -7, -6, -17, 7, -2,
                -14, -6, -4, -3, -12, -9, -9, -2, -4, -5, -17, 11, -4, -8, -5, -4,
                -13, -11, -1, -10, -5, -5, -17, 9, -10, -6, -9, -4, 2, -8, -10,
                -13, -8, -17, 8, -2, -4, -14, -13, -6, -7, -4, -5, -17, 6, 0, -5,
                -7, -6, -1, -5, -3, -17, 7, -13, -6, -3, -3, -6, -4, -17, 13, -5,
                -15, -10, -14, -11, -17, 10, -7, -6, -9, -7, -17, 7, -8, -6, -5,
                -17, 6, 0, -5, -17, 6, -5, -17, -5, -17, 1]
            )
        ),
        validate=False,
    ),
    "PAM70": SubstitutionMatrix(
        "ARNDCQEGHILKMFPSTWYVBZX*",
        _vector_to_matrix(
            np.array(
                [5, -4, -2, -1, -4, -2, -1, 0, -4, -2, -4, -4, -3, -6, 0, 1, 1, -9,
                -5, -1, -1, -1, -2, -11, 8, -3, -6, -5, 0, -5, -6, 0, -3, -6, 2,
                -2, -7, -2, -1, -4, 0, -7, -5, -4, -2, -3, -11, 6, 3, -7, -1, 0,
                -1, 1, -3, -5, 0, -5, -6, -3, 1, 0, -6, -3, -5, 5, -1, -2, -11, 6,
                -9, 0, 3, -1, -1, -5, -8, -2, -7, -10, -4, -1, -2, -10, -7, -5, 5,
                2, -3, -11, 9, -9, -9, -6, -5, -4, -10, -9, -9, -8, -5, -1, -5,
                -11, -2, -4, -8, -9, -6, -11, 7, 2, -4, 2, -5, -3, -1, -2, -9, -1,
                -3, -3, -8, -8, -4, -1, 5, -2, -11, 6, -2, -2, -4, -6, -2, -4, -9,
                -3, -2, -3, -11, -6, -4, 2, 5, -3, -11, 6, -6, -6, -7, -5, -6, -7,
                -3, 0, -3, -10, -9, -3, -1, -3, -3, -11, 8, -6, -4, -3, -6, -4, -2,
                -3, -4, -5, -1, -4, 0, 1, -3, -11, 7, 1, -4, 1, 0, -5, -4, -1, -9,
                -4, 3, -4, -4, -3, -11, 6, -5, 2, -1, -5, -6, -4, -4, -4, 0, -6,
                -4, -4, -11, 6, 0, -9, -4, -2, -1, -7, -7, -6, -1, -2, -3, -11, 10,
                -2, -5, -3, -2, -8, -7, 0, -6, -3, -3, -11, 8, -7, -4, -6, -2, 4,
                -5, -7, -9, -5, -11, 7, 0, -2, -9, -9, -3, -4, -2, -3, -11, 5, 2,
                -3, -5, -3, 0, -2, -1, -11, 6, -8, -4, -1, -1, -3, -2, -11, 13, -3,
                -10, -7, -10, -7, -11, 9, -5, -4, -7, -5, -11, 6, -5, -4, -2, -11,
                5, 1, -2, -11, 5, -3, -11, -3, -11, 1]
            )
        ),
        validate=False,
    ),
    "PAM250": SubstitutionMatrix(
        "ARNDCQEGHILKMFPSTWYVBZX*",
        _vector_to_matrix(
            np.array(
                [2, -2, 0, 0, -2, 0, 0, 1, -1, -1, -2, -1, -1, -3, 1, 1, 1, -6, -3,
                0, 0, 0, 0, -8, 6, 0, -1, -4, 1, -1, -3, 2, -2, -3, 3, 0, -4, 0, 0,
                -1, 2, -4, -2, -1, 0, -1, -8, 2, 2, -4, 1, 1, 0, 2, -2, -3, 1, -2,
                -3, 0, 1, 0, -4, -2, -2, 2, 1, 0, -8, 4, -5, 2, 3, 1, 1, -2, -4, 0,
                -3, -6, -1, 0, 0, -7, -4, -2, 3, 3, -1, -8, 12, -5, -5, -3, -3, -2,
                -6, -5, -5, -4, -3, 0, -2, -8, 0, -2, -4, -5, -3, -8, 4, 2, -1, 3,
                -2, -2, 1, -1, -5, 0, -1, -1, -5, -4, -2, 1, 3, -1, -8, 4, 0, 1,
                -2, -3, 0, -2, -5, -1, 0, 0, -7, -4, -2, 3, 3, -1, -8, 5, -2, -3,
                -4, -2, -3, -5, 0, 1, 0, -7, -5, -1, 0, 0, -1, -8, 6, -2, -2, 0,
                -2, -2, 0, -1, -1, -3, 0, -2, 1, 2, -1, -8, 5, 2, -2, 2, 1, -2, -1,
                0, -5, -1, 4, -2, -2, -1, -8, 6, -3, 4, 2, -3, -3, -2, -2, -1, 2,
                -3, -3, -1, -8, 5, 0, -5, -1, 0, 0, -3, -4, -2, 1, 0, -1, -8, 6, 0,
                -2, -2, -1, -4, -2, 2, -2, -2, -1, -8, 9, -5, -3, -3, 0, 7, -1, -4,
                -5, -2, -8, 6, 1, 0, -6, -5, -1, -1, 0, -1, -8, 2, 1, -2, -3, -1,
                0, 0, 0, -8, 3, -5, -3, 0, 0, -1, 0, -8, 17, 0, -6, -5, -6, -4, -8,
                10, -2, -3, -4, -2, -8, 4, -2, -2, -1, -8, 3, 2, -1, -8, 3, -1, -8,
                -1, -8, 1]
            )
        ),
        validate=False,
    ),
    # BLOcks SUbstitution Matrix (BLOSUM)
    "BLOSUM45": SubstitutionMatrix(
        "ARNDCQEGHILKMFPSTWYVBZX*",
        _vector_to_matrix(
            np.array(
                [ 5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -2,
                -2, 0, -1, -1, 0, -5, 7, 0, -1, -3, 1, 0, -2, 0, -3, -2, 3, -1, -2,
                -2, -1, -1, -2, -1, -2, -1, 0, -1, -5, 6, 2, -2, 0, 0, 0, 1, -2,
                -3, 0, -2, -2, -2, 1, 0, -4, -2, -3, 4, 0, -1, -5, 7, -3, 0, 2, -1,
                0, -4, -3, 0, -3, -4, -1, 0, -1, -4, -2, -3, 5, 1, -1, -5, 12, -3,
                -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -3, -2,
                -5, 6, 2, -2, 1, -2, -2, 1, 0, -4, -1, 0, -1, -2, -1, -3, 0, 4, -1,
                -5, 6, -2, 0, -3, -2, 1, -2, -3, 0, 0, -1, -3, -2, -3, 1, 4, -1,
                -5, 7, -2, -4, -3, -2, -2, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1,
                -5, 10, -3, -2, -1, 0, -2, -2, -1, -2, -3, 2, -3, 0, 0, -1, -5, 5,
                2, -3, 2, 0, -2, -2, -1, -2, 0, 3, -3, -3, -1, -5, 5, -3, 2, 1, -3,
                -3, -1, -2, 0, 1, -3, -2, -1, -5, 5, -1, -3, -1, -1, -1, -2, -1,
                -2, 0, 1, -1, -5, 6, 0, -2, -2, -1, -2, 0, 1, -2, -1, -1, -5, 8,
                -3, -2, -1, 1, 3, 0, -3, -3, -1, -5, 9, -1, -1, -3, -3, -3, -2, -1,
                -1, -5, 4, 2, -4, -2, -1, 0, 0, 0, -5, 5, -3, -1, 0, 0, -1, 0, -5,
                15, 3, -3, -4, -2, -2, -5, 8, -1, -2, -2, -1, -5, 5, -3, -3, -1,
                -5, 4, 2, -1, -5, 4, -1, -5, -1, -5, 1]
            )
        ),
        validate=False,
    ),
    "BLOSUM50": SubstitutionMatrix(
        "ARNDCQEGHILKMFPSTWYVBZX*",
        _vector_to_matrix(
            np.array(
                [5, -2, -1, -2, -1, -1, -1, 0, -2, -1, -2, -1, -1, -3, -1, 1, 0, -3,
                -2, 0, -2, -1, -1, -5, 7, -1, -2, -4, 1, 0, -3, 0, -4, -3, 3, -2,
                -3, -3, -1, -1, -3, -1, -3, -1, 0, -1, -5, 7, 2, -2, 0, 0, 0, 1,
                -3, -4, 0, -2, -4, -2, 1, 0, -4, -2, -3, 4, 0, -1, -5, 8, -4, 0, 2,
                -1, -1, -4, -4, -1, -4, -5, -1, 0, -1, -5, -3, -4, 5, 1, -1, -5,
                13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3,
                -3, -2, -5, 7, 2, -2, 1, -3, -2, 2, 0, -4, -1, 0, -1, -1, -1, -3,
                0, 4, -1, -5, 6, -3, 0, -4, -3, 1, -2, -3, -1, -1, -1, -3, -2, -3,
                1, 5, -1, -5, 8, -2, -4, -4, -2, -3, -4, -2, 0, -2, -3, -3, -4, -1,
                -2, -2, -5, 10, -4, -3, 0, -1, -1, -2, -1, -2, -3, 2, -4, 0, 0, -1,
                -5, 5, 2, -3, 2, 0, -3, -3, -1, -3, -1, 4, -4, -3, -1, -5, 5, -3,
                3, 1, -4, -3, -1, -2, -1, 1, -4, -3, -1, -5, 6, -2, -4, -1, 0, -1,
                -3, -2, -3, 0, 1, -1, -5, 7, 0, -3, -2, -1, -1, 0, 1, -3, -1, -1,
                -5, 8, -4, -3, -2, 1, 4, -1, -4, -4, -2, -5, 10, -1, -1, -4, -3,
                -3, -2, -1, -2, -5, 5, 2, -4, -2, -2, 0, 0, -1, -5, 5, -3, -2, 0,
                0, -1, 0, -5, 15, 2, -3, -5, -2, -3, -5, 8, -1, -3, -2, -1, -5, 5,
                -4, -3, -1, -5, 5, 2, -1, -5, 5, -1, -5, -1, -5, 1]
            )
        ),
        validate=False,
    ),
    "BLOSUM62": SubstitutionMatrix(
        "ARNDCQEGHILKMFPSTWYVBZX*",
        _vector_to_matrix(
            np.array(
                [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3,
                -2, 0, -2, -1, 0, -4, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3,
                -2, -1, -1, -3, -2, -3, -1, 0, -1, -4, 6, 1, -3, 0, 0, 0, 1, -3,
                -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4, 6, -3, 0, 2, -1,
                -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4, 9, -3,
                -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2,
                -4, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1,
                -4, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1,
                -4, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1,
                -4, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4, 4,
                2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4, 4, -2, 2, 0,
                -3, -2, -1, -2, -1, 1, -4, -3, -1, -4, 5, -1, -3, -1, 0, -1, -3,
                -2, -2, 0, 1, -1, -4, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4,
                6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4, 7, -1, -1, -4, -3, -2, -2,
                -1, -2, -4, 4, 1, -3, -2, -2, 0, 0, 0, -4, 5, -2, -2, 0, -1, -1, 0,
                -4, 11, 2, -3, -4, -3, -2, -4, 7, -1, -3, -2, -1, -4, 4, -3, -2,
                -1, -4, 4, 1, -1, -4, 4, -1, -4, -1, -4, 1]
            )
        ),
        validate=False,
    ),
    "BLOSUM80": SubstitutionMatrix(
        "ARNDCQEGHILKMFPSTWYVBZX*",
        _vector_to_matrix(
            np.array(
                [7, -3, -3, -3, -1, -2, -2, 0, -3, -3, -3, -1, -2, -4, -1, 2, 0, -5,
                -4, -1, -3, -2, -1, -8, 9, -1, -3, -6, 1, -1, -4, 0, -5, -4, 3, -3,
                -5, -3, -2, -2, -5, -4, -4, -2, 0, -2, -8, 9, 2, -5, 0, -1, -1, 1,
                -6, -6, 0, -4, -6, -4, 1, 0, -7, -4, -5, 5, -1, -2, -8, 10, -7, -1,
                2, -3, -2, -7, -7, -2, -6, -6, -3, -1, -2, -8, -6, -6, 6, 1, -3,
                -8, 13, -5, -7, -6, -7, -2, -3, -6, -3, -4, -6, -2, -2, -5, -5, -2,
                -6, -7, -4, -8, 9, 3, -4, 1, -5, -4, 2, -1, -5, -3, -1, -1, -4, -3,
                -4, -1, 5, -2, -8, 8, -4, 0, -6, -6, 1, -4, -6, -2, -1, -2, -6, -5,
                -4, 1, 6, -2, -8, 9, -4, -7, -7, -3, -5, -6, -5, -1, -3, -6, -6,
                -6, -2, -4, -3, -8, 12, -6, -5, -1, -4, -2, -4, -2, -3, -4, 3, -5,
                -1, 0, -2, -8, 7, 2, -5, 2, -1, -5, -4, -2, -5, -3, 4, -6, -6, -2,
                -8, 6, -4, 3, 0, -5, -4, -3, -4, -2, 1, -7, -5, -2, -8, 8, -3, -5,
                -2, -1, -1, -6, -4, -4, -1, 1, -2, -8, 9, 0, -4, -3, -1, -3, -3, 1,
                -5, -3, -2, -8, 10, -6, -4, -4, 0, 4, -2, -6, -6, -3, -8, 12, -2,
                -3, -7, -6, -4, -4, -2, -3, -8, 7, 2, -6, -3, -3, 0, -1, -1, -8, 8,
                -5, -3, 0, -1, -2, -1, -8, 16, 3, -5, -8, -5, -5, -8, 11, -3, -5,
                -4, -3, -8, 7, -6, -4, -2, -8, 6, 0, -3, -8, 6, -1, -8, -2, -8, 1]
            )
        ),
        validate=False,
    ),
    "BLOSUM90": SubstitutionMatrix(
        "ARNDCQEGHILKMFPSTWYVBZX*",
        _vector_to_matrix(
            np.array(
                [ 5, -2, -2, -3, -1, -1, -1, 0, -2, -2, -2, -1, -2, -3, -1, 1, 0, -4,
                -3, -1, -2, -1, -1, -6, 6, -1, -3, -5, 1, -1, -3, 0, -4, -3, 2, -2,
                -4, -3, -1, -2, -4, -3, -3, -2, 0, -2, -6, 7, 1, -4, 0, -1, -1, 0,
                -4, -4, 0, -3, -4, -3, 0, 0, -5, -3, -4, 4, -1, -2, -6, 7, -5, -1,
                1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5, 4, 0, -2,
                -6, 9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2,
                -4, -5, -3, -6, 7, 2, -3, 1, -4, -3, 1, 0, -4, -2, -1, -1, -3, -3,
                -3, -1, 4, -1, -6, 6, -3, -1, -4, -4, 0, -3, -5, -2, -1, -1, -5,
                -4, -3, 0, 4, -2, -6, 6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4,
                -5, -5, -2, -3, -2, -6, 8, -4, -4, -1, -3, -2, -3, -2, -2, -3, 1,
                -4, -1, 0, -2, -6, 5, 1, -4, 1, -1, -4, -3, -1, -4, -2, 3, -5, -4,
                -2, -6, 5, -3, 2, 0, -4, -3, -2, -3, -2, 0, -5, -4, -2, -6, 6, -2,
                -4, -2, -1, -1, -5, -3, -3, -1, 1, -1, -6, 7, -1, -3, -2, -1, -2,
                -2, 0, -4, -2, -1, -6, 7, -4, -3, -3, 0, 3, -2, -4, -4, -2, -6, 8,
                -2, -2, -5, -4, -3, -3, -2, -2, -6, 5, 1, -4, -3, -2, 0, -1, -1,
                -6, 6, -4, -2, -1, -1, -1, -1, -6, 11, 2, -3, -6, -4, -3, -6, 8,
                -3, -4, -3, -2, -6, 5, -4, -3, -2, -6, 4, 0, -2, -6, 4, -1, -6, -2,
                -6, 1]
            )
        ),
        validate=False,
    ),
}
# fmt: on
