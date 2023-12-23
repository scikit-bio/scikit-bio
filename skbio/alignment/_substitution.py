# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.util._decorator import experimental, classproperty, classonlymethod
# from skbio._base import SkbioObject
from skbio.stats.distance import DissimilarityMatrix


@experimental(as_of='0.5.10')
class SubstitutionMatrix(DissimilarityMatrix):
    """Matrix for scoring substitutions between characters in biological
    sequences.

    Parameters
    ----------
    alphabet : iterable
        Characters that constitute the alphabet.
    matrix : 2D array-like
        Scores of substitutions from one character (row, or axis=0) to another
        character (column, or axis=1).

    See Also
    --------
    DissimilarityMatrix

    Notes
    -----
    This class can be generalized to any combination of scalars.

    Only square matrices (i.e., numbers of rows and columns are equal) are supported

    Examples
    --------
    >>> from skbio import SubstitutionMatrix

    """

    @property
    @experimental(as_of='0.5.10')
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
        return self.ids

    @experimental(as_of='0.5.10')
    def __init__(self, data, alphabet):
        """Initialize a substitution matrix object with a score matrix and an
        alphabet.

        """
        super().__init__(data, alphabet, validate=True)

    @experimental(as_of='0.5.10')
    def to_dict(self):
        """Create a 2D dictionary from the substitution matrix.

        Returns
        -------
        dict of dict
            2D dictionary constructed from the substitution matrix.

        """
        return {id_: dict(zip(self._ids, row)) for id_, row in zip(
            self._ids, self._data)}

    @classonlymethod
    @experimental(as_of='0.5.10')
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

        Example
        -------
        >>> from skbio.alignment import SubstitutionMatrix
        >>> d = {'a': {'a': 1, 'b': 0, 'c': 0},
                 'b': {'a': 0, 'b': 1, 'c': 0},
                 'c': {'a': 0, 'b': 0, 'c': 1}}
        >>> mat = SubstitutionMatrix.from_dict(d)
        >>> mat.alphabet
        (a, b, c)
        >>> mat.data
        array([[1., 0., 0.],
               [0., 1., 0.],
               [0., 0., 1.]])
        """
        alphabet = tuple(dictionary.keys())
        alphabet_set = set(alphabet)
        idmap = {x: i for i, x in enumerate(alphabet)}
        data = np.zeros((n := len(alphabet), n))
        for i, row in enumerate(dictionary.values()):
            if set(row.keys()) != alphabet_set:
                raise ValueError('The outer and inner layers of the dictionary'
                                 ' must have the same set of keys.')
            for key, value in row.items():
                if isinstance(value, int):
                    value = float(value)
                elif not isinstance(value, float):
                    raise ValueError('Scores must be integers or floating-'
                                     'point numbers.')
                data[i][idmap[key]] = value

        return cls(data, alphabet)

    @classonlymethod
    @experimental(as_of='0.5.10')
    def identity(cls, match, mismatch, alphabet):
        """Create a substitution matrix where all matches and mismatches have
        the identical scores, respectively, regardless of the character.

        Parameters
        ----------
        match : int or float
            Score assigned to all matches.
        mismatch : int or float
            Score assigned to all mismatches.
        alphabet : iterable
            Characters that constitute the alphabet.

        Returns
        -------
        SubstitutionMatrix
            Substitution matrix constructed from the scores and alphabet.

        Example
        -------
        >>> from skbio.alignment import SubstitutionMatrix
        >>> mat = SubstitutionMatrix.identity(1, -2, 'ACGT')
        >>> mat.alphabet
        (A, C, G, T)
        >>> mat.data
        array([[1., -2., -2., -2.],
               [-2., 1., -2., -2.],
               [-2., -2., 1., -2.],
               [-2., -2., -2., 1.]])
        """
        alphabet = tuple(alphabet)
        data = np.identity(len(alphabet)) * (match - mismatch) + mismatch
        return cls(data, alphabet)


# @experimental(as_of='0.5.10')
# class SubstitutionMatrix(SkbioObject):
#     """Matrix for scoring substitutions between characters in biological
#     sequences.

#     Parameters
#     ----------
#     alphabet : iterable
#         Characters that constitute the alphabet.
#     matrix : 2D array-like
#         Scores of substitutions from one character (row, or axis=0) to another
#         character (column, or axis=1).

#     See Also
#     --------
#     DissimilarityMatrix

#     Notes
#     -----
#     This class can be generalized to any combination of scalars.

#     Only square matrices (i.e., numbers of rows and columns are equal) are supported

#     Examples
#     --------
#     >>> from skbio import SubstitutionMatrix

#     """
#     @property
#     @experimental(as_of='0.5.10')
#     def alphabet(self):
#         """Alphabet of the substitution matrix.

#         Each element (character) corresponds to one row/column in the matrix.

#         Returns
#         -------
#         list
#             Alphabet of the substitution matrix.

#         """
#         return self._alphabet

#     @property
#     @experimental(as_of='0.5.10')
#     def matrix(self):
#         """Matrix of substitution scores.

#         Each value corresponds to the score of substituting the row character
#         with the column character.

#         Returns
#         -------
#         2D np.ndarray
#             Matrix of substitution scores.

#         """
#         return self._matrix

#     @experimental(as_of='0.5.10')
#     def __init__(self, alphabet, matrix):
#         """Initialize a substitution matrix object with an alphabet and a score
#         matrix.

#         """
#         alphabet, matrix = self._validate(alphabet, matrix)
#         self._alphabet = alphabet
#         self._matrix = matrix
#         self._chardict = {c: i for i, c in enumerate(alphabet)}

#     @experimental(as_of='0.5.10')
#     def __str__(self):
#         """Return a string representation of the substitution matrix.

#         Returns
#         -------
#         str
#             String representation of the substitution matrix.

#         """
#         return ''

#     def _validate_alphabet(self, alphabet):

#         if not isinstance(alphabet, list):
#             try:
#                 alphabet = list(alphabet)
#             except TypeError:
#                 raise ValueError('Alphabet must be iterable, such as a string '
#                                  'or a list of characters.')

#         if not (num_chars := len(alphabet)):
#             raise ValueError('Alphabet must not be empty.')

#         if num_chars != len(set(alphabet)):
#             raise ValueError('Alphabet must not contain duplicated elements.')

#         return alphabet

#     def _validate_matrix(self, matrix):
#         return matrix

#     def _validate(self, alphabet, matrix):
#         alphabet = self._validate_alphabet(alphabet)
#         matrix = self._validate_matrix(matrix)
#         if len(alphabet) != matrix.shape[0]:
#             raise ValueError('Alphabet and matrix must have the same number '
#                              'of characters.')
#         return alphabet, matrix
