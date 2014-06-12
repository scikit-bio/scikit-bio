#!/usr/bin/env python
r"""
Dissimilarity and distance matrices (:mod:`skbio.core.distance`)
================================================================

.. currentmodule:: skbio.core.distance

This module provides functionality for serializing, deserializing, and
manipulating dissimilarity and distance matrices in memory. There are two
matrix classes available, `DissimilarityMatrix` and `DistanceMatrix`.
Both classes can store measures of difference/distinction between objects. A
dissimilarity/distance matrix includes both a matrix of
dissimilarities/distances (floats) between objects, as well as unique IDs
(object labels; strings) identifying each object in the matrix.

`DissimilarityMatrix` can be used to store measures of dissimilarity between
objects, and does not require that the dissimilarities are symmetric (e.g.,
dissimilarities obtained using the *Gain in PD* measure [1]_).
`DissimilarityMatrix` is a more general container to store differences than
`DistanceMatrix`.

`DistanceMatrix` has the additional requirement that the differences it
stores are symmetric (e.g., Euclidean or Hamming distances).

.. note:: `DissimilarityMatrix` can be used to store distances, but it is
   recommended to use `DistanceMatrix` to store this type of data as it
   provides an additional check for symmetry. A distance matrix is a
   dissimilarity matrix; this is modeled in the class design by having
   `DistanceMatrix` as a subclass of `DissimilarityMatrix`.

Classes
-------

.. autosummary::
   :toctree: generated/

   DissimilarityMatrix
   DistanceMatrix

Functions
---------

.. autosummary::
   :toctree: generated/

   randdm

References
----------
.. [1] Faith, D. P. (1992). "Conservation evaluation and phylogenetic
   diversity".

Examples
--------
Assume we have the following delimited text file storing distances between
three objects with IDs ``a``, ``b``, and ``c``::

    \ta\tb\tc
    a\t0.0\t0.5\t1.0
    b\t0.5\t0.0\t0.75
    c\t1.0\t0.75\t0.0

Load a distance matrix from the file:

>>> from StringIO import StringIO
>>> from skbio.core.distance import DistanceMatrix
>>> dm_f = StringIO("\ta\tb\tc\n"
...                 "a\t0.0\t0.5\t1.0\n"
...                 "b\t0.5\t0.0\t0.75\n"
...                 "c\t1.0\t0.75\t0.0\n")
>>> dm = DistanceMatrix.from_file(dm_f)
>>> print(dm)
3x3 distance matrix
IDs:
a, b, c
Data:
[[ 0.    0.5   1.  ]
 [ 0.5   0.    0.75]
 [ 1.    0.75  0.  ]]

Access the distance (scalar) between objects ``'a'`` and ``'c'``:

>>> dm['a', 'c']
1.0

Get a row vector of distances between object ``'b'`` and all other objects:

>>> dm['b']
array([ 0.5 ,  0.  ,  0.75])

numpy indexing/slicing also works as expected. Extract the third column:

>>> dm[:, 2]
array([ 1.  ,  0.75,  0.  ])

Serialize the distance matrix to delimited text file:

>>> out_f = StringIO()
>>> dm.to_file(out_f)
>>> out_f.getvalue()
'\ta\tb\tc\na\t0.0\t0.5\t1.0\nb\t0.5\t0.0\t0.75\nc\t1.0\t0.75\t0.0\n'
>>> out_f.getvalue() == dm_f.getvalue()
True

A distance matrix object can also be created from an existing ``numpy.array``
(or an array-like object, such as a nested Python list):

>>> import numpy as np
>>> data = np.array([[0.0, 0.5, 1.0],
...                  [0.5, 0.0, 0.75],
...                  [1.0, 0.75, 0.0]])
>>> ids = ["a", "b", "c"]
>>> dm_from_np = DistanceMatrix(data, ids)
>>> print(dm_from_np)
3x3 distance matrix
IDs:
a, b, c
Data:
[[ 0.    0.5   1.  ]
 [ 0.5   0.    0.75]
 [ 1.    0.75  0.  ]]
>>> dm_from_np == dm
True

IDs may be omitted when constructing a dissimilarity/distance matrix.
Monotonically-increasing integers (cast as strings) will be automatically used:

>>> dm = DistanceMatrix(data)
>>> dm.ids
('0', '1', '2')

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip
from future.utils.six import string_types

from copy import deepcopy

import numpy as np
from scipy.spatial.distance import squareform

from skbio.core.exception import (DissimilarityMatrixError,
                                  DissimilarityMatrixFormatError,
                                  DistanceMatrixError, MissingIDError)
from skbio.util.io import open_file


class DissimilarityMatrix(object):
    """Store dissimilarities between objects.

    A `DissimilarityMatrix` instance stores a square, hollow, two-dimensional
    matrix of dissimilarities between objects. Objects could be, for example,
    samples or DNA sequences. A sequence of IDs accompanies the
    dissimilarities.

    Methods are provided to load and save dissimilarity matrices from/to disk,
    as well as perform common operations such as extracting dissimilarities
    based on object ID.

    Parameters
    ----------
    data : array_like or DissimilarityMatrix
        Square, hollow, two-dimensional ``numpy.ndarray`` of dissimilarities
        (floats), or a structure that can be converted to a ``numpy.ndarray``
        using ``numpy.asarray``. Can instead be a `DissimilarityMatrix` (or
        subclass) instance, in which case the instance's data will be used.
        Data will be converted to a float ``dtype`` if necessary. A copy will
        *not* be made if already a ``numpy.ndarray`` with a float ``dtype``.
    ids : sequence of str, optional
        Sequence of strings to be used as object IDs. Must match the number of
        rows/cols in `data`. If ``None`` (the default), IDs will be
        monotonically-increasing integers cast as strings, with numbering
        starting from zero, e.g., ``('0', '1', '2', '3', ...)``.

    Attributes
    ----------
    data
    ids
    dtype
    shape
    size
    T

    See Also
    --------
    DistanceMatrix

    Notes
    -----
    The dissimilarities are stored in redundant (square-form) format [1]_.

    The data are not checked for symmetry, nor guaranteed/assumed to be
    symmetric.

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    """

    # Used in __str__
    _matrix_element_name = 'dissimilarity'

    @classmethod
    def from_file(cls, dm_f, delimiter='\t'):
        """Load dissimilarity matrix from a delimited text file or file path.

        Creates a `DissimilarityMatrix` instance from a serialized
        dissimilarity matrix stored as delimited text.

        `dm_f` can be a file-like or a file path object containing delimited
        text. The first line (header) must contain the IDs of each object. The
        subsequent lines must contain an ID followed by each dissimilarity
        (float) between the current object and all other objects, where the
        order of objects is determined by the header line.  For example, a 2x2
        dissimilarity matrix with IDs ``'a'`` and ``'b'`` might look like::

            <del>a<del>b
            a<del>0.0<del>1.0
            b<del>1.0<del>0.0

        where ``<del>`` is the delimiter between elements.

        Parameters
        ----------
        dm_f : iterable of str or str
            Iterable of strings (e.g., open file handle, file-like object, list
            of strings, etc.) or a file path (a string) containing a serialized
            dissimilarity matrix.
        delimiter : str, optional
            String delimiting elements in `dm_f`.

        Returns
        -------
        DissimilarityMatrix
            Instance of type `cls` containing the parsed contents of `dm_f`.

        Notes
        -----
        Whitespace-only lines can occur anywhere throughout the "file" and are
        ignored. Lines starting with ``#`` are treated as comments and ignored.
        These comments can only occur *before* the ID header.

        IDs will have any leading/trailing whitespace removed when they are
        parsed.

        .. note::
            File-like objects passed to this method will not be closed upon the
            completion of the parsing, it is responsibility of the owner of the
            object to perform this operation.

        """
        # We aren't using np.loadtxt because it uses *way* too much memory
        # (e.g, a 2GB matrix eats up 10GB, which then isn't freed after parsing
        # has finished). See:
        # http://mail.scipy.org/pipermail/numpy-tickets/2012-August/006749.html

        with open_file(dm_f, 'U') as dm_f:

            # We use iter() as we want to take a single pass over the
            # iterable and maintain our current position after finding
            # the header (mainly necessary for something like a list
            # of strings).
            dm_f = iter(dm_f)

            # Strategy:
            #   - find the header
            #   - initialize an empty ndarray
            #   - for each row of data in the input file:
            #     - populate the corresponding row in the ndarray with floats

            ids = cls._parse_ids(dm_f, delimiter)
            num_ids = len(ids)
            data = np.empty((num_ids, num_ids), dtype=np.float64)

            # curr_row_idx keeps track of the row index within the data matrix.
            # We're not using enumerate() because there may be
            # empty/whitespace-only lines throughout the data matrix. We want
            # to ignore those and only count the actual rows of data.
            curr_row_idx = 0
            for line in dm_f:
                line = line.strip()

                if not line:
                    continue
                elif curr_row_idx >= num_ids:
                    # We've hit a nonempty line after we already filled the
                    # data matrix. Raise an error because we shouldn't ignore
                    # extra data.
                    raise DissimilarityMatrixFormatError(
                        "Encountered extra rows without corresponding IDs in"
                        " the header.")

                tokens = line.split(delimiter)

                # -1 because the first element contains the current ID.
                if len(tokens) - 1 != num_ids:
                    raise DissimilarityMatrixFormatError(
                        "There are %d values in row number %d, which is not"
                        " equal to the number of IDs in the header (%d)."
                        % (len(tokens) - 1, curr_row_idx + 1, num_ids))

                curr_id = tokens[0].strip()
                expected_id = ids[curr_row_idx]
                if curr_id == expected_id:
                    data[curr_row_idx, :] = np.asarray(tokens[1:], dtype=float)
                else:
                    raise DissimilarityMatrixFormatError(
                        "Encountered mismatched IDs while parsing the "
                        "dissimilarity matrix file. Found '%s' but expected "
                        "'%s'. Please ensure that the IDs match between the "
                        "dissimilarity matrix header (first row) and the row "
                        "labels (first column)." % (curr_id, expected_id))

                curr_row_idx += 1

        if curr_row_idx != num_ids:
            raise DissimilarityMatrixFormatError(
                "Expected %d row(s) of data, but found %d." % (num_ids,
                                                               curr_row_idx))

        return cls(data, ids)

    def __init__(self, data, ids=None):
        if isinstance(data, DissimilarityMatrix):
            data = data.data
        data = np.asarray(data, dtype='float')

        if ids is None:
            ids = (str(i) for i in range(data.shape[0]))
        ids = tuple(ids)

        self._validate(data, ids)

        self._data = data
        self._ids = ids
        self._id_index = self._index_list(self._ids)

    @property
    def data(self):
        """Array of dissimilarities.

        A square, hollow, two-dimensional ``numpy.ndarray`` of dissimilarities
        (floats). A copy is *not* returned.

        Notes
        -----
        This property is not writeable.

        """
        return self._data

    @property
    def ids(self):
        """Tuple of object IDs.

        A tuple of strings, one for each object in the dissimilarity matrix.

        Notes
        -----
        This property is writeable, but the number of new IDs must match the
        number of objects in `data`.

        """
        return self._ids

    @ids.setter
    def ids(self, ids_):
        ids_ = tuple(ids_)
        self._validate(self.data, ids_)
        self._ids = ids_
        self._id_index = self._index_list(self._ids)

    @property
    def dtype(self):
        """Data type of the dissimilarities."""
        return self.data.dtype

    @property
    def shape(self):
        """Two-element tuple containing the dissimilarity matrix dimensions.

        Notes
        -----
        As the dissimilarity matrix is guaranteed to be square, both tuple
        entries will always be equal.

        """
        return self.data.shape

    @property
    def size(self):
        """Total number of elements in the dissimilarity matrix.

        Notes
        -----
        Equivalent to ``self.shape[0] * self.shape[1]``.

        """
        return self.data.size

    @property
    def T(self):
        """Transpose of the dissimilarity matrix.

        See Also
        --------
        transpose

        """
        return self.transpose()

    def transpose(self):
        """Return the transpose of the dissimilarity matrix.

        Notes
        -----
        A deep copy is returned.

        Returns
        -------
        DissimilarityMatrix
            Transpose of the dissimilarity matrix. Will be the same type as
            `self`.

        """
        return self.__class__(self.data.T.copy(), deepcopy(self.ids))

    def redundant_form(self):
        """Return an array of dissimilarities in redundant format.

        As this is the native format that the dissimilarities are stored in,
        this is simply an alias for `data`.

        Returns
        -------
        ndarray
            Two-dimensional ``numpy.ndarray`` of dissimilarities in redundant
            format.

        Notes
        -----
        Redundant format is described in [1]_.

        Does *not* return a copy of the data.

        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

        """
        return self.data

    def copy(self):
        """Return a deep copy of the dissimilarity matrix.

        Returns
        -------
        DissimilarityMatrix
            Deep copy of the dissimilarity matrix. Will be the same type as
            `self`.

        """
        # We deepcopy IDs in case the tuple contains mutable objects at some
        # point in the future.
        return self.__class__(self.data.copy(), deepcopy(self.ids))

    def __str__(self):
        """Return a string representation of the dissimilarity matrix.

        Summary includes matrix dimensions, a (truncated) list of IDs, and
        (truncated) array of dissimilarities.

        Returns
        -------
        str
            String representation of the dissimilarity matrix.

        .. shownumpydoc

        """
        return '%dx%d %s matrix\nIDs:\n%s\nData:\n' % (
            self.shape[0], self.shape[1], self._matrix_element_name,
            self._pprint_ids()) + str(self.data)

    def __eq__(self, other):
        """Compare this dissimilarity matrix to another for equality.

        Two dissimilarity matrices are equal if they have the same shape, IDs
        (in the same order!), and have data arrays that are equal.

        Checks are *not* performed to ensure that `other` is a
        `DissimilarityMatrix` instance.

        Parameters
        ----------
        other : DissimilarityMatrix
            Dissimilarity matrix to compare to for equality.

        Returns
        -------
        bool
            ``True`` if `self` is equal to `other`, ``False`` otherwise.

        .. shownumpydoc

        """
        equal = True

        # The order these checks are performed in is important to be as
        # efficient as possible. The check for shape equality is not strictly
        # necessary as it should be taken care of in np.array_equal, but I'd
        # rather explicitly bail before comparing IDs or data. Use array_equal
        # instead of (a == b).all() because of this issue:
        #     http://stackoverflow.com/a/10582030
        try:
            if self.shape != other.shape:
                equal = False
            elif self.ids != other.ids:
                equal = False
            elif not np.array_equal(self.data, other.data):
                equal = False
        except AttributeError:
            equal = False

        return equal

    def __ne__(self, other):
        """Determine whether two dissimilarity matrices are not equal.

        Parameters
        ----------
        other : DissimilarityMatrix
            Dissimilarity matrix to compare to.

        Returns
        -------
        bool
            ``True`` if `self` is not equal to `other`, ``False`` otherwise.

        See Also
        --------
        __eq__

        .. shownumpydoc

        """
        return not self == other

    def __getitem__(self, index):
        """Slice into dissimilarity data by object ID or numpy indexing.

        Extracts data from the dissimilarity matrix by object ID, a pair of
        IDs, or numpy indexing/slicing.

        Parameters
        ----------
        index : str, two-tuple of str, or numpy index
            `index` can be one of the following forms: an ID, a pair of IDs, or
            a numpy index.

            If `index` is a string, it is assumed to be an ID and a
            ``numpy.ndarray`` row vector is returned for the corresponding ID.
            Note that the ID's row of dissimilarities is returned, *not* its
            column. If the matrix is symmetric, the two will be identical, but
            this makes a difference if the matrix is asymmetric.

            If `index` is a two-tuple of strings, each string is assumed to be
            an ID and the corresponding matrix element is returned that
            represents the dissimilarity between the two IDs. Note that the
            order of lookup by ID pair matters if the matrix is asymmetric: the
            first ID will be used to look up the row, and the second ID will be
            used to look up the column. Thus, ``dm['a', 'b']`` may not be the
            same as ``dm['b', 'a']`` if the matrix is asymmetric.

            Otherwise, `index` will be passed through to
            ``DissimilarityMatrix.data.__getitem__``, allowing for standard
            indexing of a ``numpy.ndarray`` (e.g., slicing).

        Returns
        -------
        ndarray or scalar
            Indexed data, where return type depends on the form of `index` (see
            description of `index` for more details).

        Raises
        ------
        MissingIDError
            If the ID(s) specified in `index` are not in the dissimilarity
            matrix.

        Notes
        -----
        The lookup based on ID(s) is quick.

        .. shownumpydoc

        """
        if isinstance(index, string_types):
            if index in self._id_index:
                return self.data[self._id_index[index]]
            else:
                raise MissingIDError(index)
        elif self._is_id_pair(index):
            for id_ in index:
                if id_ not in self._id_index:
                    raise MissingIDError(id_)
            return self.data[self._id_index[index[0]],
                             self._id_index[index[1]]]
        else:
            return self.data.__getitem__(index)

    def to_file(self, out_f, delimiter='\t'):
        """Save the dissimilarity matrix to file in delimited text format.

        Parameters
        ----------
        out_f : file-like object or filename
            File-like object to write serialized data to, or name of
            file. If it's a file-like object, it must have a ``write``
            method, and it won't be closed. Else, it is opened and
            closed after writing.
        delimiter : str, optional
            Delimiter used to separate elements in output format.

        See Also
        --------
        from_file

        """
        with open_file(out_f, 'w') as out_f:
            formatted_ids = self._format_ids(delimiter)
            out_f.write(formatted_ids)
            out_f.write('\n')

            for id_, vals in zip(self.ids, self.data):
                out_f.write(id_)
                out_f.write(delimiter)
                out_f.write(delimiter.join(np.asarray(vals, dtype=np.str)))
                out_f.write('\n')

    @staticmethod
    def _parse_ids(dm_f, delimiter):
        header_line = None

        for line in dm_f:
            line = line.strip()

            if line and not line.startswith('#'):
                header_line = line
                break

        if header_line is None:
            raise DissimilarityMatrixFormatError(
                "Could not find a header line containing IDs in the "
                "dissimilarity matrix file. Please verify that the file is "
                "not empty.")
        else:
            return [e.strip() for e in header_line.split(delimiter)]

    def _validate(self, data, ids):
        """Validate the data array and IDs.

        Checks that the data is at least 1x1 in size, 2D, square, hollow, and
        contains only floats. Also checks that IDs are unique and that the
        number of IDs matches the number of rows/cols in the data array.

        Subclasses can override this method to perform different/more specific
        validation (e.g., see `DistanceMatrix`).

        Notes
        -----
        Accepts arguments instead of inspecting instance attributes to avoid
        creating an invalid dissimilarity matrix before raising an error.
        Otherwise, the invalid dissimilarity matrix could be used after the
        exception is caught and handled.

        """
        num_ids = len(ids)

        if 0 in data.shape:
            raise DissimilarityMatrixError("Data must be at least 1x1 in "
                                           "size.")
        elif len(data.shape) != 2:
            raise DissimilarityMatrixError("Data must have exactly two "
                                           "dimensions.")
        elif data.shape[0] != data.shape[1]:
            raise DissimilarityMatrixError("Data must be square (i.e., have "
                                           "the same number of rows and "
                                           "columns).")
        elif data.dtype != np.double:
            raise DissimilarityMatrixError("Data must contain only floating "
                                           "point values.")
        elif np.trace(data) != 0:
            raise DissimilarityMatrixError("Data must be hollow (i.e., the "
                                           "diagonal can only contain zeros).")
        elif num_ids != len(set(ids)):
            raise DissimilarityMatrixError("IDs must be unique.")
        elif num_ids != data.shape[0]:
            raise DissimilarityMatrixError("The number of IDs must match the "
                                           "number of rows/columns in the "
                                           "data.")

    def _index_list(self, list_):
        return {id_: idx for idx, id_ in enumerate(list_)}

    def _is_id_pair(self, index):
        return (isinstance(index, tuple) and
                len(index) == 2 and
                all(map(lambda e: isinstance(e, string_types), index)))

    def _format_ids(self, delimiter):
        return delimiter.join([''] + list(self.ids))

    def _pprint_ids(self, max_chars=80, delimiter=', ', suffix='...',):
        # Adapted from http://stackoverflow.com/a/250373
        ids_str = delimiter.join(self.ids)

        if len(ids_str) > max_chars:
            truncated = ids_str[:max_chars + 1].split(delimiter)[0:-1]
            ids_str = delimiter.join(truncated) + delimiter + suffix

        return ids_str


class DistanceMatrix(DissimilarityMatrix):
    """Store distances between objects.

    A `DistanceMatrix` is a `DissimilarityMatrix` with the additional
    requirement that the matrix data is symmetric. There are additional methods
    made available that take advantage of this symmetry.

    See Also
    --------
    DissimilarityMatrix

    Notes
    -----
    The distances are stored in redundant (square-form) format [1]_. To
    facilitate use with other scientific Python routines (e.g., scipy), the
    distances can be retrieved in condensed (vector-form) format using
    `condensed_form`.

    `DistanceMatrix` only requires that the distances it stores are symmetric.
    Checks are *not* performed to ensure the other three metric properties
    hold (non-negativity, identity of indiscernibles, and triangle inequality)
    [2]_. Thus, a `DistanceMatrix` instance can store distances that are not
    metric.

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html
    .. [2] http://planetmath.org/metricspace

    """

    # Override here, used in superclass __str__
    _matrix_element_name = 'distance'

    def condensed_form(self):
        """Return an array of distances in condensed format.

        Returns
        -------
        ndarray
            One-dimensional ``numpy.ndarray`` of distances in condensed format.

        Notes
        -----
        Condensed format is described in [1]_.

        The conversion is not a constant-time operation, though it should be
        relatively quick to perform.

        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

        """
        return squareform(self._data, force='tovector', checks=False)

    def permute(self, condensed=False):
        """Randomly permute both rows and columns in the matrix.

        Randomly permutes the ordering of rows and columns in the matrix. The
        same permutation is applied to both rows and columns in order to
        maintain symmetry and hollowness. Only the rows/columns in the distance
        matrix are permuted; the IDs are *not* permuted.

        Parameters
        ----------
        condensed : bool, optional
            If ``True``, return the permuted distance matrix in condensed
            format. Otherwise, return the permuted distance matrix as a new
            ``DistanceMatrix`` instance.

        Returns
        -------
        DistanceMatrix or ndarray
            Permuted distances as a new ``DistanceMatrix`` or as a ``ndarray``
            in condensed format.

        See Also
        --------
        condensed_form

        Notes
        -----
        This method does not modify the distance matrix that it is called on.
        It is more efficient to pass ``condensed=True`` than permuting the
        distance matrix and then converting to condensed format.

        """
        order = np.random.permutation(self.shape[0])
        permuted = self._data[order][:, order]

        if condensed:
            return squareform(permuted, force='tovector', checks=False)
        else:
            return self.__class__(permuted, self.ids)

    def _validate(self, data, ids):
        """Validate the data array and IDs.

        Overrides the superclass `_validate`. Performs a check for symmetry in
        addition to the checks performed in the superclass.

        """
        super(DistanceMatrix, self)._validate(data, ids)

        if (data.T != data).any():
            raise DistanceMatrixError("Data must be symmetric.")


def randdm(num_objects, ids=None, constructor=None, random_fn=None):
    """Generate a distance matrix populated with random distances.

    Using the default `random_fn`, distances are randomly drawn from a uniform
    distribution over ``[0, 1)``.

    Regardless of `random_fn`, the resulting distance matrix is guaranteed to
    be symmetric and hollow.

    Parameters
    ----------
    num_objects : int
        The number of objects in the resulting distance matrix. For example, if
        `num_objects` is 3, a 3x3 distance matrix will be returned.
    ids : sequence of str or None, optional
        A sequence of strings to be used as IDs. ``len(ids)`` must be equal to
        `num_objects`. If not provided, IDs will be monotonically-increasing
        integers cast as strings (numbering starts at 1). For example,
        ``('1', '2', '3')``.
    constructor : type, optional
        `DissimilarityMatrix` or subclass constructor to use when creating the
        random distance matrix. The returned distance matrix will be of this
        type. If ``None`` (the default), a `DistanceMatrix` instance will be
        returned.
    random_fn : function, optional
        Function to generate random values. `random_fn` must accept two
        arguments (number of rows and number of columns) and return a 2D
        ``numpy.ndarray`` of floats (or something that can be cast to float).
        If ``None`` (the default), ``numpy.random.rand`` will be used.

    Returns
    -------
    DissimilarityMatrix
        `DissimilarityMatrix` (or subclass) instance of random distances. Type
        depends on `constructor`.

    See Also
    --------
    numpy.random.rand

    """
    if constructor is None:
        constructor = DistanceMatrix
    if random_fn is None:
        random_fn = np.random.rand

    data = np.tril(random_fn(num_objects, num_objects), -1)
    data += data.T

    if not ids:
        ids = map(str, range(1, num_objects + 1))

    return constructor(data, ids)
