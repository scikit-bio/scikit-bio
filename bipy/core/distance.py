#!/usr/bin/env python
r"""
Distance matrices (:mod:`bipy.core.distance`)
=============================================

.. currentmodule:: bipy.core.distance

This module provides functionality for serializing, deserializing, and
manipulating distance matrices in memory. There are multiple distance matrix
classes available, where the appropriate class to use depends on the nature of
the distances you wish to store.

A distance matrix includes both a matrix of distances (floats) between objects,
as well as IDs (labels) identifying each object in the matrix.

Classes
-------

.. autosummary::
   :toctree: generated/

   DistanceMatrix
   SymmetricDistanceMatrix

Functions
---------

.. autosummary::
   :toctree: generated/

   random_distance_matrix

Examples
--------
Assume we have the following delimited text file::

    \ta\tb\tc
    a\t0.0\t0.5\t1.0
    b\t0.5\t0.0\t0.75
    c\t1.0\t0.75\t0.0

Load a distance matrix from the file:

>>> from StringIO import StringIO
>>> from bipy.core.distance import DistanceMatrix
>>> dm_f = StringIO("\ta\tb\tc\n"
...                 "a\t0.0\t0.5\t1.0\n"
...                 "b\t0.5\t0.0\t0.75\n"
...                 "c\t1.0\t0.75\t0.0\n")
>>> dm = DistanceMatrix.from_file(dm_f)
>>> print dm
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

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from copy import deepcopy
from itertools import izip

import numpy as np
from scipy.spatial.distance import squareform


class DistanceMatrixError(Exception):
    r"""General error for distance matrix validation failures."""
    pass


class MissingIDError(Exception):
    r"""Error for ID lookup that doesn't exist in the distance matrix."""

    def __init__(self, id_):
        super(MissingIDError, self).__init__()
        self.args = ("The ID '%s' is not in the distance matrix." % id_,)


class DistanceMatrixFormatError(Exception):
    r"""Error for reporting issues in distance matrix file format.

    Typically used during parsing.

    """
    pass


class IDMismatchError(Exception):
    r"""Error for reporting a mismatch between IDs.

    Typically used during parsing.

    """

    def __init__(self, actual, expected):
        super(IDMismatchError, self).__init__()
        self.args = ("Encountered mismatched IDs while parsing the distance "
                     "matrix file. Found '%s' but expected '%s'. Please "
                     "ensure that the IDs match between the distance matrix "
                     "header (first row) and the row labels (first column)." %
                     (actual, expected),)


class MissingHeaderError(Exception):
    r"""Error for reporting a missing ID header line during parsing."""

    def __init__(self):
        super(MissingHeaderError, self).__init__()
        self.args = ("Could not find a header line containing IDs in the "
                     "distance matrix file. Please verify that the file is "
                     "not empty.",)


class MissingDataError(Exception):
    r"""Error for reporting missing data lines during parsing."""

    def __init__(self, actual, expected):
        super(MissingDataError, self).__init__()
        self.args = ("Expected %d row(s) of data, but found %d." % (expected,
                                                                    actual),)


class DistanceMatrix(object):
    """Store distances between objects and object IDs.

    A `DistanceMatrix` instance stores a square, hollow, two-dimensional matrix
    of distances between objects. Objects could be, for example, samples or DNA
    sequences. A sequence of IDs accompanies the distances.

    Methods are provided to load and save distance matrices from/to disk, as
    well as perform common operations such as extracting distances based on ID.

    Parameters
    ----------
    data : array_like or DistanceMatrix
        Square, hollow, two-dimensional ``numpy.ndarray`` of distances
        (floats), or a structure that can be converted to a ``numpy.ndarray``
        using ``numpy.asarray``. Can instead be a `DistanceMatrix` instance, in
        which case the distance matrix's data will be used. Data will be
        converted to a float ``dtype`` if necessary. A copy will *not* be made
        if already a ``numpy.ndarray`` with a float ``dtype``.
    ids : sequence of str
        Sequence of strings to be used as object IDs. Must match the number of
        rows/cols in `data`.

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
    SymmetricDistanceMatrix

    Notes
    -----
    The distances are stored in redundant (square-form) format [1]_.

    The data are not checked for symmetry, nor guaranteed/assumed to be
    symmetric.

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    """

    @classmethod
    def from_file(cls, dm_f, delimiter='\t'):
        r"""Load distance matrix from delimited text file.

        Creates a `DistanceMatrix` instance from a serialized distance
        matrix stored as delimited text.

        `dm_f` must be a file-like object containing delimited text. The first
        line (header) must contain the IDs of each object. The subsequent lines
        must contain an ID followed by each distance (float) between the
        current object and all other objects, where the order of objects is
        determined by the header line.

        For example, a 2x2 distance matrix with IDs ``'a'`` and ``'b'`` might
        look like::

            <tab>a<tab>b
            a<tab>0.0<tab>1.0
            b<tab>1.0<tab>0.0

        where ``<tab>`` is the delimiter between elements.

        Parameters
        ----------
        dm_f : iterable of str
            Iterable of strings (e.g., open file handle, file-like object, list
            of strings, etc.) containing a serialized distance matrix.
        delimiter : str, optional
            String delimiting elements in `dm_f`.

        Returns
        -------
        DistanceMatrix
            Instance of type `cls` containing the parsed contents of `dm_f`.

        Notes
        -----
        Whitespace-only lines can occur anywhere throughout the "file" and are
        ignored. Lines starting with ``#`` are treated as comments and ignored.
        These comments can only occur *before* the ID header.

        IDs will have any leading/trailing whitespace removed when they are
        parsed.

        """
        # We aren't using np.loadtxt because it uses *way* too much memory
        # (e.g, a 2GB matrix eats up 10GB, which then isn't freed after parsing
        # has finished). See:
        # http://mail.scipy.org/pipermail/numpy-tickets/2012-August/006749.html

        # Strategy:
        #     - find the header
        #     - initialize an empty ndarray
        #     - for each row of data in the input file:
        #         - populate the corresponding row in the ndarray with floats

        # We use iter() as we want to take a single pass over the iterable and
        # maintain our current position after finding the header (mainly
        # necessary for something like a list of strings).
        dm_f = iter(dm_f)
        ids = cls._parse_ids(dm_f, delimiter)
        num_ids = len(ids)
        data = np.empty((num_ids, num_ids), dtype='float')

        # curr_row_idx keeps track of the row index within the data matrix.
        # We're not using enumerate() because there may be
        # empty/whitespace-only lines throughout the data matrix. We want to
        # ignore those and only count the actual rows of data.
        curr_row_idx = 0
        for line in dm_f:
            line = line.strip()

            if not line:
                continue
            elif curr_row_idx >= num_ids:
                # We've hit a nonempty line after we already filled the data
                # matrix. Raise an error because we shouldn't ignore extra
                # data.
                raise DistanceMatrixFormatError(
                    "Encountered extra rows without corresponding IDs in the "
                    "header.")

            tokens = line.split(delimiter)

            # -1 because the first element contains the current ID.
            if len(tokens) - 1 != num_ids:
                raise DistanceMatrixFormatError(
                    "There are %d values in row number %d, which is not equal "
                    "to the number of IDs in the header (%d)."
                    % (len(tokens) - 1, curr_row_idx + 1, num_ids))

            curr_id = tokens[0].strip()
            expected_id = ids[curr_row_idx]
            if curr_id == expected_id:
                data[curr_row_idx, :] = np.asarray(tokens[1:], dtype='float')
            else:
                raise IDMismatchError(curr_id, expected_id)

            curr_row_idx += 1

        if curr_row_idx != num_ids:
            raise MissingDataError(curr_row_idx, num_ids)

        return cls(data, ids)

    def __init__(self, data, ids):
        if isinstance(data, DistanceMatrix):
            data = data.data
        data = np.asarray(data, dtype='float')

        ids = tuple(ids)
        self._validate(data, ids)

        self._data = data
        self._ids = ids
        self._id_index = self._index_list(self._ids)

    @property
    def data(self):
        r"""Array of distances.

        A square, hollow, two-dimensional ``numpy.ndarray`` of distances
        (floats). A copy is *not* returned.

        Notes
        -----
        This property is not writeable.

        """
        return self._data

    @property
    def ids(self):
        r"""Tuple of object IDs.

        A tuple of strings, one for each object in the distance matrix.

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
        r"""Data type of the distances."""
        return self.data.dtype

    @property
    def shape(self):
        r"""Two-element tuple containing the distance matrix dimensions.

        Notes
        -----
        As the distance matrix is guaranteed to be square, both tuple entries
        will be equal.

        """
        return self.data.shape

    @property
    def size(self):
        r"""Total number of elements in the distance matrix.

        Notes
        -----
        Equivalent to ``self.shape[0] * self.shape[1]``.

        """
        return self.data.size

    @property
    def T(self):
        r"""Transpose of the distance matrix.

        See Also
        --------
        transpose

        """
        return self.transpose()

    def transpose(self):
        r"""Return the transpose of the distance matrix.

        Notes
        -----
        A deep copy is returned.

        Returns
        -------
        DistanceMatrix
            Transpose of the distance matrix. Will be the same type as `self`.

        """
        return self.__class__(self.data.T.copy(), deepcopy(self.ids))

    def redundant_form(self):
        r"""Return an array of distances in redundant format.

        As this is the native format that the distances are stored in, this is
        simply an alias for `data`.

        Returns
        -------
        ndarray
            Two-dimensional ``numpy.ndarray`` of distances in redundant format.

        Notes
        -----
        Redundant and condensed formats are described in [1]_.

        Does *not* return a copy of the data.

        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

        """
        return self.data

    def copy(self):
        r"""Return a deep copy of the distance matrix.

        Returns
        -------
        DistanceMatrix
            Deep copy of the distance matrix. Will be the same type as `self`.

        """
        # We deepcopy IDs in case the tuple contains mutable objects at some
        # point in the future.
        return self.__class__(self.data.copy(), deepcopy(self.ids))

    def __str__(self):
        r"""Return a string representation of the distance matrix.

        Summary includes matrix dimensions, a (truncated) list of IDs, and
        (truncated) array of distances.

        Returns
        -------
        str
            String representation (summary) of the distance matrix.

        """
        return '%dx%d distance matrix\nIDs:\n%s\nData:\n' % (
            self.shape[0], self.shape[1],
            self._pprint_ids()) + str(self.data)

    def __eq__(self, other):
        r"""Compare this distance matrix to another for equality.

        Two distance matrices are equal if they have the same shape, IDs (in
        the same order!), and have data arrays that are equal.

        Checks are *not* performed to ensure that `other` is a `DistanceMatrix`
        instance.

        Parameters
        ----------
        other : DistanceMatrix
            Distance matrix to compare to for equality.

        Returns
        -------
        bool
            ``True`` if `self` is equal to `other`, ``False`` otherwise.

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
        r"""Determine whether two distance matrices are not equal.

        Parameters
        ----------
        other : DistanceMatrix
            Distance matrix to compare to.

        Returns
        -------
        bool
            ``True`` if `self` is not equal to `other`, ``False`` otherwise.

        See Also
        --------
        __eq__

        """
        return not self == other

    def __getitem__(self, index):
        r"""Slice into distance data by ID or numpy indexing.

        Extracts data from the distance matrix by ID, a pair of IDs, or numpy
        indexing/slicing.

        Parameters
        ----------
        index : str, two-tuple of str, or numpy index
            `index` can be one of the following forms: an ID, a pair of IDs, or
            a numpy index.

            If `index` is a string, it is assumed to be an ID and a
            ``numpy.ndarray`` row vector is returned for the corresponding ID.
            Note that the ID's row of distances is returned, *not* its column.
            If the matrix is symmetric, the two will be identical, but this
            makes a difference if the matrix is asymmetric.

            If `index` is a two-tuple of strings, each string is assumed to be
            an ID and the corresponding matrix element is returned that
            represents the distance between the two IDs. Note that the order of
            lookup by ID pair matters if the matrix is asymmetric: the first ID
            will be used to look up the row, and the second ID will be used to
            look up the column. Thus, ``dm['a', 'b']`` may not be the same as
            ``dm['b', 'a']`` if the matrix is asymmetric.

            Otherwise, `index` will be passed through to
            ``DistanceMatrix.data.__getitem__``, allowing for standard indexing
            of a ``numpy.ndarray`` (e.g., slicing).

        Returns
        -------
        ndarray or scalar
            Indexed data, where return type depends on the form of `index` (see
            description of `index` for more details).

        Raises
        ------
        MissingIDError
            If the ID(s) specified in `index` are not in the distance matrix.

        Notes
        -----
        The lookup based on ID(s) is quick.

        """
        if isinstance(index, basestring):
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
        r"""Save the distance matrix to file in delimited text format.

        See Also
        --------
        from_file

        Parameters
        ----------
        out_f : file-like object
            File-like object to write serialized data to. Must have a ``write``
            method. It is the caller's responsibility to close `out_f` when
            done (if necessary).
        delimiter : str, optional
            Delimiter used to separate elements in output format.

        """
        formatted_ids = self._format_ids(delimiter)
        out_f.write(formatted_ids)
        out_f.write('\n')

        for id_, vals in izip(self.ids, self.data):
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
            raise MissingHeaderError
        else:
            return map(lambda e: e.strip(), header_line.split(delimiter))

    def _validate(self, data, ids):
        r"""Validate the data array and IDs.

        Checks that the data is at least 1x1 in size, 2D, square, hollow, and
        contains only floats. Also checks that IDs are unique and that the
        number of IDs matches the number of rows/cols in the data array.

        Subclasses can override this method to perform different/more specific
        validation (e.g., see `SymmetricDistanceMatrix`).

        Notes
        -----
        Accepts arguments instead of inspecting instance attributes to avoid
        creating an invalid distance matrix before raising an error. Otherwise,
        the invalid distance matrix could be used after the exception is
        caught and handled.

        """
        num_ids = len(ids)

        if 0 in data.shape:
            raise DistanceMatrixError("Data must be at least 1x1 in size.")
        elif len(data.shape) != 2:
            raise DistanceMatrixError("Data must have exactly two dimensions.")
        elif data.shape[0] != data.shape[1]:
            raise DistanceMatrixError("Data must be square (i.e., have the "
                                      "same number of rows and columns).")
        elif data.dtype != np.double:
            raise DistanceMatrixError("Data must contain only floating point "
                                      "values.")
        elif np.trace(data) != 0:
            raise DistanceMatrixError("Data must be hollow (i.e., the "
                                      "diagonal can only contain zeros).")
        elif num_ids != len(set(ids)):
            raise DistanceMatrixError("IDs must be unique.")
        elif num_ids != data.shape[0]:
            raise DistanceMatrixError("The number of IDs must match the "
                                      "number of rows/columns in the data.")

    def _index_list(self, list_):
        return {id_: idx for idx, id_ in enumerate(list_)}

    def _is_id_pair(self, index):
        return (isinstance(index, tuple) and
                len(index) == 2 and
                all(map(lambda e: isinstance(e, basestring), index)))

    def _format_ids(self, delimiter):
        return delimiter.join([''] + list(self.ids))

    def _pprint_ids(self, max_chars=80, delimiter=', ', suffix='...',):
        # Adapted from http://stackoverflow.com/a/250373
        ids_str = delimiter.join(self.ids)

        if len(ids_str) > max_chars:
            truncated = ids_str[:max_chars + 1].split(delimiter)[0:-1]
            ids_str = delimiter.join(truncated) + delimiter + suffix

        return ids_str


class SymmetricDistanceMatrix(DistanceMatrix):
    r"""Store symmetric distance data.

    A `SymmetricDistanceMatrix` is a `DistanceMatrix` with the additional
    requirement that the matrix data is symmetric. There are additional methods
    made available that take advantage of this symmetry.

    See Also
    --------
    DistanceMatrix

    Notes
    -----
    The distances are stored in redundant (square-form) format [1]_. To
    facilitate use with other scientific Python routines (e.g., scipy), the
    distances can be retrieved in condensed (vector-form) format using
    `condensed_form`.

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    """

    def condensed_form(self):
        r"""Return an array of distances in condensed format.

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
        return squareform(self.data, force='tovector')

    def _validate(self, data, ids):
        r"""Validate the data array and IDs.

        Overrides the superclass `_validate`. Performs a check for symmetry in
        addition to the checks performed in the superclass.

        """
        super(SymmetricDistanceMatrix, self)._validate(data, ids)

        if (data.T != data).any():
            raise DistanceMatrixError("Data must be symmetric.")


def random_distance_matrix(num_objects, ids=None, constructor=DistanceMatrix,
                           random_fn=np.random.rand):
    r"""Generate a distance matrix populated with random distances.

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
        `DistanceMatrix` or subclass constructor to use when creating the
        random distance matrix. The returned distance matrix will be of this
        type.
    random_fn : function, optional
        Function to generate random values. `random_fn` must accept two
        arguments (number of rows and number of columns) and return a 2D
        ``numpy.ndarray`` of floats (or something that can be cast to float).

    Returns
    -------
    DistanceMatrix
        `DistanceMatrix` (or subclass) instance of random distances. Type
        depends on `constructor`.

    See Also
    --------
    numpy.random.rand

    """
    data = np.tril(random_fn(num_objects, num_objects), -1)
    data += data.T

    if not ids:
        ids = map(str, range(1, num_objects + 1))

    return constructor(data, ids)
