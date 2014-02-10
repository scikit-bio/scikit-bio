#!/usr/bin/env python
"""Core distance matrix API."""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from copy import deepcopy
from itertools import izip

import numpy as np
from scipy.spatial.distance import is_valid_dm, squareform


class DistanceMatrixError(Exception):
    """General error for distance matrix validation failures."""
    pass


class MissingSampleIDError(Exception):
    """Error for sample ID lookup that doesn't exist in the distance matrix."""
    pass


class DistanceMatrixFormatError(Exception):
    """Error for reporting issues in distance matrix file format.

    Typically used during parsing.

    """
    pass


class SampleIDMismatchError(Exception):
    """Error for reporting a mismatch between sample IDs.

    Typically used during parsing.

    """

    def __init__(self, actual, expected):
        super(SampleIDMismatchError, self).__init__()
        self.args = ("Encountered mismatched sample IDs while parsing the "
                     "distance matrix file. Found '%s' but expected '%s'. "
                     "Please ensure that the sample IDs match between the "
                     "distance matrix header (first row) and the row labels "
                     "(first column)." % (actual, expected),)


class MissingHeaderError(Exception):
    """Error for reporting a missing sample ID header line during parsing."""

    def __init__(self):
        super(MissingHeaderError, self).__init__()
        self.args = ("Could not find a header line containing sample IDs in "
                     "the distance matrix file. Please verify that the file "
                     "is not empty.",)


class MissingDataError(Exception):
    """Error for reporting missing data lines during parsing."""

    def __init__(self, actual, expected):
        super(MissingDataError, self).__init__()
        self.args = ("Expected %d row(s) of data, but found %d." % (expected,
                                                                    actual),)


class DistanceMatrix(object):
    """Encapsulate a 2D array of distances (floats) and sample IDs (labels).

    A ``DistanceMatrix`` instance contains a square, hollow 2D numpy
    ``ndarray`` of distances (floats) between samples. A tuple of sample IDs
    (typically strings) accompanies the raw distance data. Methods are provided
    to load and save distance matrices, as well as perform common operations
    such as extracting distances based on sample ID.

    The ``ndarray`` of distances can be accessed via the ``data`` property. The
    tuple of sample IDs can be accessed via the ``sample_ids`` property.
    ``sample_ids`` is also writeable, though the new sample IDs must match the
    number of samples in ``data``.

    The distances are stored in redundant (square-form) format. For more
    details on redundant format, see:
    http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    Note that the data is not checked for symmetry, nor guaranteed/assumed to
    be symmetric. See the subclass ``SymmetricDistanceMatrix`` for validation
    and operations that take advantage of matrix symmetry.

    """

    @classmethod
    def from_file(cls, dm_f, delimiter='\t'):
        """Load distance matrix from delimited text file.

        Given an open file-like object ``dm_f``, an instance of this class is
        returned based on the parsed file contents.

        The file must contain delimited text (controlled via ``delimiter``).
        The first line must contain all sample IDs, where each ID is separated
        by ``delimiter``. The subsequent lines must contain a sample ID
        followed by each distance (float) between the sample and all other
        samples.

        For example, a 2x2 distance matrix with samples ``'a'`` and ``'b'``
        would look like:

        <tab>a<tab>b
        a<tab>0.0<tab>1.0
        b<tab>1.0<tab>0.0

        where ``<tab>`` is the delimiter between elements.

        Whitespace-only lines can occur anywhere throughout the file and are
        ignored. Lines starting with # are treated as comments and ignored.
        These comments can only occur *before* the sample ID header.

        Sample IDs will have any leading/trailing whitespace removed.

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
        sids, header_idx = cls._parse_sample_ids(dm_f, delimiter)

        # If we have a list of strings, we need to "seek" to the start of the
        # data section. If we have a file-like object, we should already be
        # there.
        try:
            data_f = dm_f[(header_idx + 1):]
        except AttributeError:
            data_f = dm_f

        num_sids = len(sids)
        data = np.empty((num_sids, num_sids), dtype='float')

        # curr_row_idx keeps track of the row index within the data matrix.
        # We're not using enumerate() because there may be
        # empty/whitespace-only lines throughout the data matrix. We want to
        # ignore those and only count the actual rows of data.
        curr_row_idx = 0
        for line in data_f:
            line = line.strip()

            if not line:
                continue
            elif curr_row_idx >= num_sids:
                # We've hit a nonempty line after we already filled the data
                # matrix. Raise an error because we shouldn't ignore extra
                # data.
                raise DistanceMatrixFormatError(
                    "Encountered extra rows without corresponding sample IDs "
                    "in the header.")

            tokens = line.split(delimiter)

            # -1 because the first column contains the sample ID.
            if len(tokens) - 1 != num_sids:
                raise DistanceMatrixFormatError(
                    "There are %d values in row number %d, which is not equal "
                    "to the number of sample IDs in the header (%d)."
                    % (len(tokens) - 1, curr_row_idx + 1, num_sids))

            curr_sid = tokens[0].strip()
            expected_sid = sids[curr_row_idx]
            if curr_sid == expected_sid:
                data[curr_row_idx, :] = np.asarray(tokens[1:], dtype='float')
            else:
                raise SampleIDMismatchError(curr_sid, expected_sid)

            curr_row_idx += 1

        if curr_row_idx != num_sids:
            raise MissingDataError(curr_row_idx, num_sids)

        return cls(data, sids)

    def __init__(self, data, sample_ids):
        """Construct a ``DistanceMatrix`` instance.

        Arguments:
        data -- a square, hollow 2D ``numpy.ndarray`` of distances (floats), or
            a structure that can be converted to a ``numpy.ndarray`` using
            ``numpy.asarray``. Data will be converted to a float ``dtype`` if
            necessary. A copy will *not* be made if already an ``ndarray`` with
            a float ``dtype``
        sample_ids -- a sequence of strings to be used as sample labels. Must
            match the number of rows/cols in ``data``

        """
        data = np.asarray(data, dtype='float')
        sample_ids = tuple(sample_ids)
        self._validate(data, sample_ids)

        self._data = data
        self._sample_ids = sample_ids
        self._sample_index = self._index_list(self._sample_ids)

    @property
    def data(self):
        """Return a ``numpy.ndarray`` of distances.

        A copy is *not* returned. This property is not writeable.

        """
        return self._data

    @property
    def sample_ids(self):
        """Return a tuple of sample IDs.

        This property is writeable, but the number of new sample IDs must match
        the number of samples in ``data``.

        """
        return self._sample_ids

    @sample_ids.setter
    def sample_ids(self, sample_ids_):
        sample_ids_ = tuple(sample_ids_)
        self._validate(self.data, sample_ids_)
        self._sample_ids = sample_ids_
        self._sample_index = self._index_list(self._sample_ids)

    @property
    def dtype(self):
        """Return the ``dtype`` of the underlying ``numpy.ndarray``."""
        return self.data.dtype

    @property
    def shape(self):
        """Return a two-element tuple containing the array dimensions.

        As the distance matrix is guaranteed to be square, both tuple entries
        will be equal.

        """
        return self.data.shape

    @property
    def num_samples(self):
        """Returns the number of samples (i.e. number of rows or columns)."""
        return len(self.sample_ids)

    @property
    def size(self):
        """Return the total number of elements in the distance matrix.

        Equivalent to ``self.shape[0] * self.shape[1]``.

        """
        return self.data.size

    @property
    def T(self):
        """Return the transpose of the distance matrix."""
        return self.transpose()

    def transpose(self):
        """Return the transpose of the distance matrix."""
        return self.__class__(self.data.T.copy(), deepcopy(self.sample_ids))

    def redundant_form(self):
        """Return a 2D ``numpy.ndarray`` of distances.

        As this is the native format that the distances are stored in, this is
        simply an alias for ``self.data``.

        Does *not* return a copy of the data.

        For more details on redundant and condensed formats, see:
        http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

        """
        return self.data

    def copy(self):
        """Return a deep copy of the distance matrix."""
        # We deepcopy sample IDs in case the tuple contains mutable objects at
        # some point in the future.
        return self.__class__(self.data.copy(), deepcopy(self.sample_ids))

    def __str__(self):
        """Return a string representation of the distance matrix.

        Summary includes matrix dimensions, a (truncated) list of sample IDs,
        and (truncated) array of distances.

        """
        return '%dx%d distance matrix\nSample IDs:\n%s\nData:\n' % (
            self.shape[0], self.shape[1],
            self._pprint_sample_ids()) + str(self.data)

    def __eq__(self, other):
        """Return ``True`` if this distance matrix is equal to the other.

        Two distance matrices are equal if they have the same shape, sample IDs
        (in the same order!), and have data arrays that are equal.

        Checks are *not* performed to ensure that ``other`` is a
        ``DistanceMatrix`` instance.

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
            elif self.sample_ids != other.sample_ids:
                equal = False
            elif not np.array_equal(self.data, other.data):
                equal = False
        except AttributeError:
            equal = False

        return equal

    def __ne__(self, other):
        """Return ``True`` if this distance matrix and the other are not equal.

        See ``__eq__`` for more details.

        """
        return not self == other

    def __getitem__(self, index):
        """Slice into the data by sample ID or numpy indexing.

        If ``index`` is a string, it is assumed to be a sample ID and a
        ``numpy.ndarray`` row vector is returned for the corresponding sample.
        Note that the sample's row of distances is returned, not its column. If
        the matrix is symmetric, the two will be identical, but this makes a
        difference if the matrix is asymmetric.

        The lookup based on sample ID is quick. ``MissingSampleIDError`` is
        raised if the sample does not exist.

        Otherwise, ``index`` will be passed through to
        ``DistanceMatrix.data.__getitem__``, allowing for standard indexing of
        a numpy ``ndarray`` (e.g., slicing).

        Arguments:
        index -- the sample ID whose row of distances will be returned, or the
            index to be passed through to the underlying data matrix

        """
        if isinstance(index, basestring):
            if index in self._sample_index:
                return self.data[self._sample_index[index]]
            else:
                raise MissingSampleIDError("The sample ID '%s' is not in the "
                                           "distance matrix." % index)
        else:
            return self.data.__getitem__(index)

    def to_file(self, out_f, delimiter='\t'):
        """Save the distance matrix to file in delimited text format.

        See ``from_file`` for more details on the file format.

        Arguments:
        out_f -- file-like object to write to
        delimiter -- delimiter used to separate elements in output format

        """
        formatted_sids = self._format_sample_ids(delimiter)
        out_f.write(formatted_sids)
        out_f.write('\n')

        for sid, vals in izip(self.sample_ids, self.data):
            out_f.write(sid)
            out_f.write(delimiter)
            out_f.write(delimiter.join(np.asarray(vals, dtype=np.str)))
            out_f.write('\n')

    @staticmethod
    def _parse_sample_ids(dm_f, delimiter):
        header_line = None

        for idx, line in enumerate(dm_f):
            line = line.strip()

            if line and not line.startswith('#'):
                header_line = line
                break

        if header_line is None:
            raise MissingHeaderError
        else:
            return map(lambda e: e.strip(), header_line.split(delimiter)), idx

    def _validate(self, data, sample_ids):
        """Validate the data array and sample IDs.

        Checks that the data is at least 1x1 in size, 2D, square, hollow, and
        contains only floats. Also checks that sample IDs are unique and that
        the number of sample IDs matches the number of rows/cols in the data
        array.

        Subclasses can override this method to perform different/more specific
        validation (e.g., see ``SymmetricDistanceMatrix``).

        Accepts arguments instead of inspecting instance attributes because we
        don't want to create an invalid distance matrix before raising an error
        (because then it could be used after the exception is caught).

        """
        num_sids = len(sample_ids)

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
                                      "diagonal can only contain zeros.)")
        elif num_sids != len(set(sample_ids)):
            raise DistanceMatrixError("Sample IDs must be unique.")
        elif num_sids != data.shape[0]:
            raise DistanceMatrixError("The number of sample IDs must match "
                                      "the number of rows/columns in the "
                                      "data.")

    def _index_list(self, list_):
        return {id_: idx for idx, id_ in enumerate(list_)}

    def _format_sample_ids(self, delimiter):
        return delimiter.join([''] + list(self.sample_ids))

    def _pprint_sample_ids(self, max_chars=80, delimiter=', ', suffix='...',):
        """Adapted from http://stackoverflow.com/a/250373"""
        sids_str = delimiter.join(self.sample_ids)

        if len(sids_str) > max_chars:
            truncated = sids_str[:max_chars + 1].split(delimiter)[0:-1]
            sids_str = delimiter.join(truncated) + delimiter + suffix

        return sids_str


class SymmetricDistanceMatrix(DistanceMatrix):
    """Represent a symmetric distance matrix.

    A ``SymmetricDistanceMatrix`` is a ``DistanceMatrix`` with the additional
    requirement that the matrix data is symmetric. There are additional
    operations made available that take advantage of this symmetry.

    The distances are stored in redundant (square-form) format (same as
    ``DistanceMatrix``). To facilitate use with other scientific Python
    routines (e.g., scipy), the distances can be retrieved in condensed
    (vector-form) format using ``condensed_form``. For more details on
    redundant and condensed formats, see:
    http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    """

    def condensed_form(self):
        """Return a 1D ``numpy.ndarray`` vector of distances.

        The conversion is not a constant-time operation, though it should be
        relatively quick to perform.

        For more details on redundant and condensed formats, see:
        http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

        """
        return squareform(self.data, force='tovector')

    def _validate(self, data, sample_ids):
        """Validate the data array and sample IDs.

        Overrides the superclass ``_validate``. Performs a check for symmetry
        in addition to the checks performed in the superclass.

        """
        super(SymmetricDistanceMatrix, self)._validate(data, sample_ids)

        if (data.T != data).any():
            raise DistanceMatrixError("Data must be symmetric.")


def random_distance_matrix(num_samples, sample_ids=None,
                           constructor=DistanceMatrix,
                           random_fn=np.random.rand):
    """Return a distance matrix populated with random distances.

    Using the default ``random_fn``, distances are randomly drawn from a
    uniform distribution over ``[0, 1)`` (see ``numpy.random.rand`` for more
    details).

    Regardless of the ``random_fn`` that is used to populate the matrix, the
    returned distance matrix is guaranteed to be symmetric and hollow.

    Arguments:
    num_samples -- the number of samples in the resulting distance matrix. For
        example, if ``num_samples`` is 3, a 3x3 distance matrix will be
        returned
    sample_ids -- a sequence of strings to be used as sample IDs.
        ``len(sample_ids)`` must be equal to ``num_samples``. If not provided,
        sample IDs will be monotonically-increasing integers cast as strings
        (numbering starts at 1). For example, ``('1', '2', '3')``
    constructor -- ``DistanceMatrix`` or subclass constructor to use when
        creating the distance matrix. The returned distance matrix will be of
        this type
    random_fn -- function to generate random values. Function must accept two
        arguments (number of rows and number of columns) and return a 2D
        ``numpy.ndarray`` of floats (or something that can be casted to float)

    """
    data = np.tril(random_fn(num_samples, num_samples), -1)
    data += data.T

    if not sample_ids:
        sample_ids = map(str, range(1, num_samples + 1))

    return constructor(data, sample_ids)
