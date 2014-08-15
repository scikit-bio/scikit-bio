# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import zip
from future.utils.six import StringIO, string_types

import csv
from copy import deepcopy

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

from skbio.io.util import open_file
from skbio.stats import p_value_to_str


class DissimilarityMatrixError(Exception):
    """General error for dissimilarity matrix validation failures."""
    pass


class DistanceMatrixError(DissimilarityMatrixError):
    """General error for distance matrix validation failures."""
    pass


class MissingIDError(DissimilarityMatrixError):
    """Error for ID lookup that doesn't exist in the dissimilarity matrix."""

    def __init__(self, missing_id):
        super(MissingIDError, self).__init__()
        self.args = ("The ID '%s' is not in the dissimilarity matrix." %
                     missing_id,)


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

    def index(self, lookup_id):
        """Return the index of the specified ID.

        Parameters
        ----------
        lookup_id : str
            ID whose index will be returned.

        Returns
        -------
        int
            Row/column index of `lookup_id`.

        Raises
        ------
        MissingIDError
            If `lookup_id` is not in the dissimilarity matrix.

        """
        if lookup_id in self:
            return self._id_index[lookup_id]
        else:
            raise MissingIDError(lookup_id)

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

    def filter(self, ids, strict=True):
        """Filter the dissimilarity matrix by IDs.

        Parameters
        ----------
        ids : iterable of str
            IDs to retain. May not contain duplicates or be empty. Each ID must
            be present in the dissimilarity matrix.
        strict : bool, optional
            If `strict` is ``True`` and an ID that is not found in the distance
            matrix is found in `ids`, a ``MissingIDError`` exception will be
            raised, otherwise the ID will be ignored.

        Returns
        -------
        DissimilarityMatrix
            Filtered dissimilarity matrix containing only the IDs specified in
            `ids`. IDs will be in the same order as they appear in `ids`.

        Raises
        ------
        MissingIDError
            If an ID in `ids` is not in the object's list of IDs.
        """
        if strict:
            idxs = [self.index(id_) for id_ in ids]
        else:
            # get the indices to slice the inner numpy array
            idxs = []
            # save the IDs that were found in the distance matrix
            found_ids = []
            for id_ in ids:
                try:
                    idxs.append(self.index(id_))
                    found_ids.append(id_)
                except MissingIDError:
                    pass
            ids = found_ids

        filtered_data = self._data[idxs][:, idxs]
        return self.__class__(filtered_data, ids)

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

    def __contains__(self, lookup_id):
        """Check if the specified ID is in the dissimilarity matrix.

        Parameters
        ----------
        lookup_id : str
            ID to search for.

        Returns
        -------
        bool
            ``True`` if `lookup_id` is in the dissimilarity matrix, ``False``
            otherwise.

        See Also
        --------
        index

        .. shownumpydoc

        """
        return lookup_id in self._id_index

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
            return self.data[self.index(index)]
        elif self._is_id_pair(index):
            return self.data[self.index(index[0]), self.index(index[1])]
        else:
            return self.data.__getitem__(index)

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


class CategoricalStats(object):
    """Base class for categorical statistical methods.

    Categorical statistical methods generally test for significant differences
    between discrete groups of objects, as determined by a categorical variable
    (grouping vector).

    See Also
    --------
    ANOSIM
    PERMANOVA

    """

    short_method_name = ''
    long_method_name = ''
    test_statistic_name = ''

    def __init__(self, distance_matrix, grouping, column=None):
        if not isinstance(distance_matrix, DistanceMatrix):
            raise TypeError("Input must be a DistanceMatrix.")

        if isinstance(grouping, pd.DataFrame):
            if column is None:
                raise ValueError("Must provide a column name if supplying a "
                                 "data frame.")
            else:
                grouping = self._df_to_vector(distance_matrix, grouping,
                                              column)
        elif column is not None:
            raise ValueError("Must provide a data frame if supplying a column "
                             "name.")

        if len(grouping) != distance_matrix.shape[0]:
            raise ValueError("Grouping vector size must match the number of "
                             "IDs in the distance matrix.")

        # Find the group labels and convert grouping to an integer vector
        # (factor).
        groups, grouping = np.unique(grouping, return_inverse=True)

        if len(groups) == len(grouping):
            raise ValueError("All values in the grouping vector are unique. "
                             "This method cannot operate on a grouping vector "
                             "with only unique values (e.g., there are no "
                             "'within' distances because each group of "
                             "objects contains only a single object).")
        if len(groups) == 1:
            raise ValueError("All values in the grouping vector are the same. "
                             "This method cannot operate on a grouping vector "
                             "with only a single group of objects (e.g., "
                             "there are no 'between' distances because there "
                             "is only a single group).")

        self._dm = distance_matrix
        self._grouping = grouping
        self._groups = groups
        self._tri_idxs = np.triu_indices(self._dm.shape[0], k=1)

    def _df_to_vector(self, distance_matrix, df, column):
        """Return a grouping vector from a data frame column.

        Parameters
        ----------
        distance_marix : DistanceMatrix
            Distance matrix whose IDs will be mapped to group labels.
        df : pandas.DataFrame
            Data frame (indexed by distance matrix ID).
        column : str
            Column name in `df` containing group labels.

        Returns
        -------
        list
            Grouping vector (vector of labels) based on the IDs in
            `distance_matrix`. Each ID's label is looked up in the data frame
            under the column specified by `column`.

        Raises
        ------
        ValueError
            If `column` is not in the data frame, or a distance matrix ID is
            not in the data frame.

        """
        if column not in df:
            raise ValueError("Column '%s' not in data frame." % column)

        grouping = df.loc[distance_matrix.ids, column]
        if grouping.isnull().any():
            raise ValueError("One or more IDs in the distance matrix are not "
                             "in the data frame.")
        return grouping.tolist()

    def __call__(self, permutations=999):
        """Execute the statistical method.

        Parameters
        ----------
        permutations : int, optional
            Number of permutations to use when calculating statistical
            significance. Must be >= 0. If 0, the resulting p-value will be
            ``None``.

        Returns
        -------
        CategoricalStatsResults
            Results of the method, including test statistic and p-value.

        .. shownumpydoc

        """
        if permutations < 0:
            raise ValueError("Number of permutations must be greater than or "
                             "equal to zero.")

        stat = self._run(self._grouping)

        p_value = None
        if permutations > 0:
            perm_stats = np.empty(permutations, dtype=np.float64)

            for i in range(permutations):
                perm_grouping = np.random.permutation(self._grouping)
                perm_stats[i] = self._run(perm_grouping)

            p_value = ((perm_stats >= stat).sum() + 1) / (permutations + 1)

        return CategoricalStatsResults(self.short_method_name,
                                       self.long_method_name,
                                       self.test_statistic_name,
                                       self._dm.shape[0], self._groups, stat,
                                       p_value, permutations)

    def _run(self, grouping):
        raise NotImplementedError("Subclasses must implement _run().")


class CategoricalStatsResults(object):
    """Statistical method results container.

    Stores the results of running a `CategoricalStats` method a single time,
    and provides a way to format the results.

    Attributes
    ----------
    short_method_name
    long_method_name
    test_statistic_name
    sample_size
    groups
    statistic
    p_value
    permutations

    Notes
    -----
    Users will generally not directly instantiate objects of this class. The
    various categorical statistical methods will return an object of this type
    when they are run.

    """

    def __init__(self, short_method_name, long_method_name,
                 test_statistic_name, sample_size, groups, statistic, p_value,
                 permutations):
        self.short_method_name = short_method_name
        self.long_method_name = long_method_name
        self.test_statistic_name = test_statistic_name
        self.sample_size = sample_size
        self.groups = groups
        self.statistic = statistic
        self.p_value = p_value
        self.permutations = permutations

    def __str__(self):
        """Return pretty-print (fixed width) string."""
        rows = (self._format_header(), self._format_data())

        max_widths = []
        for col_idx in range(len(rows[0])):
            max_widths.append(max(map(lambda e: len(e[col_idx]), rows)))

        results = []
        for row in rows:
            padded_row = []
            for col_idx, val in enumerate(row):
                padded_row.append(val.rjust(max_widths[col_idx]))
            results.append('  '.join(padded_row))

        return '\n'.join(results) + '\n'

    def _repr_html_(self):
        """Return a string containing an HTML table of results.

        This method will be called within the IPython Notebook instead of
        __repr__ to display results.

        """
        header = self._format_header()
        data = self._format_data()
        return pd.DataFrame([data[1:]], columns=header[1:],
                            index=[data[0]])._repr_html_()

    def summary(self, delimiter='\t'):
        """Return a formatted summary of results as a string.

        The string is formatted as delimited text.

        Parameters
        ----------
        delimiter : str, optional
            String to delimit fields by in formatted output. Default is tab
            (TSV).

        Returns
        -------
        str
            Delimited-text summary of results.

        """
        summary = StringIO()
        csv_writer = csv.writer(summary, delimiter=delimiter,
                                lineterminator='\n')
        csv_writer.writerow(self._format_header())
        csv_writer.writerow(self._format_data())
        return summary.getvalue()

    def _format_header(self):
        return ('Method name', 'Sample size', 'Number of groups',
                self.test_statistic_name, 'p-value', 'Number of permutations')

    def _format_data(self):
        p_value_str = p_value_to_str(self.p_value, self.permutations)

        return (self.short_method_name, '%d' % self.sample_size,
                '%d' % len(self.groups), str(self.statistic), p_value_str,
                '%d' % self.permutations)
