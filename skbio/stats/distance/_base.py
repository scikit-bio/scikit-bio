# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
from copy import deepcopy

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

from skbio._base import SkbioObject
from skbio.stats._misc import _pprint_strs
from skbio.util import find_duplicates
from skbio.util._decorator import classonlymethod
from skbio.util._misc import resolve_key
from skbio.util._plotting import PlottableMixin

from ._utils import is_symmetric_and_hollow
from ._utils import distmat_reorder, distmat_reorder_condensed


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
        self.args = ("The ID '%s' is not in the dissimilarity matrix." % missing_id,)


class DissimilarityMatrix(SkbioObject, PlottableMixin):
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
        using ``numpy.asarray`` or a one-dimensional vector of dissimilarities
        (floats), as defined by `scipy.spatial.distance.squareform`. Can
        instead be a `DissimilarityMatrix` (or subclass) instance,
        in which case the instance's data will be used.
        Data will be converted to a float ``dtype`` if necessary. A copy will
        *not* be made if already a ``numpy.ndarray`` with a float ``dtype``.
    ids : sequence of str, optional
        Sequence of strings to be used as object IDs. Must match the number of
        rows/cols in `data`. If ``None`` (the default), IDs will be
        monotonically-increasing integers cast as strings, with numbering
        starting from zero, e.g., ``('0', '1', '2', '3', ...)``.
    validate : bool, optional
        If `validate` is ``True`` (the default) and data is not a
        DissimilarityMatrix object, the input data will be validated.

    See Also
    --------
    DistanceMatrix
    scipy.spatial.distance.squareform

    Notes
    -----
    The dissimilarities are stored in redundant (square-form) format [1]_.

    The data are not checked for symmetry, nor guaranteed/assumed to be
    symmetric.

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    """

    default_write_format = "lsmat"
    # Used in __str__
    _matrix_element_name = "dissimilarity"

    def __init__(self, data, ids=None, validate=True):
        validate_full = validate
        validate_shape = False
        validate_ids = False

        if isinstance(data, DissimilarityMatrix):
            if isinstance(data, self.__class__):
                # Never validate when copying from an object
                # of the same type
                # We should be able to assume it is already
                # in a good state.
                validate_full = False
                validate_shape = False
                # but do validate ids, if redefining them
                validate_ids = False if ids is None else True
            ids = data.ids if ids is None else ids
            data = data.data

        # It is necessary to standardize the representation of the .data
        # attribute of this object. The input types might be list, tuple,
        # np.array, or possibly some other object type. Generally, this
        # normalization of type will require a copy of data. For example,
        # moving from a Python type representation (e.g., [[0, 1], [1, 0]])
        # requires casting all of the values to numpy types, which is handled
        # as an implicit copy via np.asarray. However, these copies are
        # unnecessary if the data object is already a numpy array. np.asarray
        # is smart enough to not copy the data, however if a dtype change is
        # requested it will. The following block of code limits the use of
        # np.asarray to situations where the data are (a) not already a numpy
        # array or (b) the data are not a single or double precision numpy
        # data type.
        _issue_copy = True
        if isinstance(data, np.ndarray):
            if data.dtype in (np.float32, np.float64):
                _issue_copy = False

        if _issue_copy:
            data = np.asarray(data, dtype="float")

        if data.ndim == 1:
            # We can assume squareform will return a symmetric square matrix
            # so no need for full validation.
            # Still do basic checks (e.g. zero length)
            # and id validation
            data = squareform(data, force="tomatrix", checks=False)
            validate_full = False
            validate_shape = True
            validate_ids = True

        if ids is None:
            ids = (str(i) for i in range(data.shape[0]))
            # I just created the ids, so no need to re-validate them
            validate_ids = False
        ids = tuple(ids)

        if validate_full:
            self._validate(data, ids)
        else:
            if validate_shape:
                self._validate_shape(data)
            if validate_ids:
                self._validate_ids(data, ids)

        self._data = data
        self._ids = ids
        self._id_index = self._index_list(self._ids)

    @classonlymethod
    def from_iterable(cls, iterable, metric, key=None, keys=None):
        """Create DissimilarityMatrix from an iterable given a metric.

        Parameters
        ----------
        iterable : iterable
            Iterable containing objects to compute pairwise dissimilarities on.
        metric : callable
            A function that takes two arguments and returns a float
            representing the dissimilarity between the two arguments.
        key : callable or metadata key, optional
            A function that takes one argument and returns a string
            representing the id of the element in the dissimilarity matrix.
            Alternatively, a key to a `metadata` property if it exists for
            each element in the `iterable`. If None, then default ids will be
            used.
        keys : iterable, optional
            An iterable of the same length as `iterable`. Each element will be
            used as the respective key.

        Returns
        -------
        DissimilarityMatrix
            The `metric` applied to all pairwise elements in the `iterable`.

        Raises
        ------
        ValueError
            If `key` and `keys` are both provided.

        """
        iterable = list(iterable)
        if key is not None and keys is not None:
            raise ValueError("Cannot use both `key` and `keys` at the same" " time.")

        keys_ = None
        if key is not None:
            keys_ = [resolve_key(e, key) for e in iterable]
        elif keys is not None:
            keys_ = keys

        dm = np.empty((len(iterable),) * 2)
        for i, a in enumerate(iterable):
            for j, b in enumerate(iterable):
                dm[i, j] = metric(a, b)

        return cls(dm, keys_)

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
        self._validate_ids(self.data, ids_)
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
        # Note: Skip validation, since we assume self was already validated
        return self.__class__(self.data.T.copy(), deepcopy(self.ids), validate=False)

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
        # Note: Skip validation, since we assume self was already validated
        return self.__class__(self.data.copy(), deepcopy(self.ids), validate=False)

    def rename(self, mapper, strict=True):
        """Rename IDs in the dissimilarity matrix.

        Parameters
        ----------
        mapper : dict or callable
            A dictionary or function that maps current IDs to new IDs.
        strict : bool, optional
           If ``True`` (default), every ID in the matrix must be included in
           ``mapper``. If ``False``, only the specified IDs will be renamed.

        Raises
        ------
        ValueError
            If ``mapper`` does not contain all of the same IDs in the matrix
            whereas in strict mode.

        Examples
        --------
        >>> from skbio import DistanceMatrix
        >>> dm = DistanceMatrix([[0, 1], [1, 0]], ids=['a', 'b'])
        >>> dm.rename({'a': 'x', 'b': 'y'})
        >>> print(dm.ids)
        ('x', 'y')

        """
        if isinstance(mapper, dict):
            if strict and not set(self.ids).issubset(mapper):
                raise ValueError(
                    "The IDs in mapper do not include all IDs in the matrix."
                )
            new_ids = [mapper.get(x, x) for x in self.ids]
        else:
            new_ids = [mapper(x) for x in self.ids]
        self.ids = new_ids

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
        if tuple(self._ids) == tuple(ids):
            return self.__class__(self._data, self._ids)

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

        # Note: Skip validation, since we assume self was already validated
        # But ids are new, so validate them explicitly
        filtered_data = distmat_reorder(self._data, idxs)
        self._validate_ids(filtered_data, ids)
        return self.__class__(filtered_data, ids, validate=False)

    def _stable_order(self, ids):
        """Obtain a stable ID order with respect to self.

        Parameters
        ----------
        ids : Iterable of ids
            The IDs to establish a stable ordering for.

        Returns
        -------
        np.array, dtype=int
            The corresponding index values

        """
        id_order = sorted(self._id_index[i] for i in ids)
        return np.array(id_order, dtype=int)

    def within(self, ids):
        """Obtain all the distances among the set of IDs.

        Parameters
        ----------
        ids : Iterable of str
            The IDs to obtain distances for. All pairs of distances are
            returned such that, if provided ['a', 'b', 'c'], the distances
            for [('a', 'a'), ('a', 'b'), ('a', 'c'), ('b', 'a'), ('b', 'b'),
            ('b', 'c'), ('c', 'a'), ('c', 'b'), ('c', 'c')] are gathered.

        Returns
        -------
        pd.DataFrame
            (i, j, value) representing the source ID ("i"), the target ID ("j")
            and the distance ("value").

        Raises
        ------
        MissingIDError
            If an ID(s) specified is not in the dissimilarity matrix.

        Notes
        -----
        Order of the return items is stable, meaning that requesting IDs
        ['a', 'b'] is equivalent to ['b', 'a']. The order is with respect
        to the order of the .ids attribute of self.

        Examples
        --------
        >>> from skbio.stats.distance import DissimilarityMatrix
        >>> dm = DissimilarityMatrix([[0, 1, 2, 3, 4], [1, 0, 1, 2, 3],
        ...                           [2, 1, 0, 1, 2], [3, 2, 1, 0, 1],
        ...                           [4, 3, 2, 1, 0]],
        ...                          ['A', 'B', 'C', 'D', 'E'])
        >>> dm.within(['A', 'B', 'C'])
           i  j  value
        0  A  A    0.0
        1  A  B    1.0
        2  A  C    2.0
        3  B  A    1.0
        4  B  B    0.0
        5  B  C    1.0
        6  C  A    2.0
        7  C  B    1.0
        8  C  C    0.0

        """
        ids = set(ids)
        not_present = ids - set(self._id_index)
        if not_present:
            raise MissingIDError(
                "At least one ID (e.g., '%s') was not " "found." % not_present.pop()
            )

        return self._subset_to_dataframe(ids, ids)

    def between(self, from_, to_, allow_overlap=False):
        """Obtain the distances between the two groups of IDs.

        Parameters
        ----------
        from_ : Iterable of str
            The IDs to obtain distances from. Distances from all pairs of IDs
            in from and to will be obtained.
        to_ : Iterable of str
            The IDs to obtain distances to. Distances from all pairs of IDs
            in to and from will be obtained.

        allow_overlap : bool, optional
            If True, allow overlap in the IDs of from and to (which would in
            effect be collecting the within distances). Default is False.

        Returns
        -------
        pd.DataFrame
            (i, j, value) representing the source ID ("i"), the target ID ("j")
            and the distance ("value").

        Raises
        ------
        MissingIDError
            If an ID(s) specified is not in the dissimilarity matrix.

        Notes
        -----
        Order of the return items is stable, meaning that requesting IDs
        ['a', 'b'] is equivalent to ['b', 'a']. The order is with respect to
        the .ids attribute of self.

        Examples
        --------
        >>> from skbio.stats.distance import DissimilarityMatrix
        >>> dm = DissimilarityMatrix([[0, 1, 2, 3, 4], [1, 0, 1, 2, 3],
        ...                           [2, 1, 0, 1, 2], [3, 2, 1, 0, 1],
        ...                           [4, 3, 2, 1, 0]],
        ...                          ['A', 'B', 'C', 'D', 'E'])
        >>> dm.between(['A', 'B'], ['C', 'D', 'E'])
           i  j  value
        0  A  C    2.0
        1  A  D    3.0
        2  A  E    4.0
        3  B  C    1.0
        4  B  D    2.0
        5  B  E    3.0

        """
        from_ = set(from_)
        to_ = set(to_)

        all_ids = from_ | to_
        not_present = all_ids - set(self._id_index)
        if not_present:
            raise MissingIDError(
                "At least one ID (e.g., '%s') was not " "found." % not_present.pop()
            )

        overlapping = from_ & to_
        if not allow_overlap and overlapping:
            raise KeyError(
                "At least one ID overlaps in from_ and to_ "
                "(e.g., '%s'). This constraint can removed with "
                "allow_overlap=True." % overlapping.pop()
            )

        return self._subset_to_dataframe(from_, to_)

    def _subset_to_dataframe(self, i_ids, j_ids):
        """Extract a subset of self and express as a DataFrame.

        Parameters
        ----------
        i_ids : Iterable of str
            The "from" IDs.
        j_ids : Iterable of str
            The "to" IDs.

        Notes
        -----
        ID membership is not tested by this private method, and it is assumed
        the caller has asserted the IDs are present.

        Returns
        -------
        pd.DataFrame
            (i, j, value) representing the source ID ("i"), the target ID ("j")
            and the distance ("value").

        """
        i_indices = self._stable_order(i_ids)
        j_indices = self._stable_order(j_ids)

        j_length = len(j_indices)
        j_labels = tuple([self.ids[j] for j in j_indices])

        i = []
        j = []

        # np.hstack([]) throws a ValueError. However, np.hstack([np.array([])])
        # is valid and returns an empty array. Accordingly, an empty array is
        # included here so that np.hstack works in the event that either i_ids
        # or j_ids is empty.
        values = [np.array([])]
        for i_idx in i_indices:
            i.extend([self.ids[i_idx]] * j_length)
            j.extend(j_labels)

            subset = self._data[i_idx, j_indices]
            values.append(subset)

        i = pd.Series(i, name="i", dtype=str)
        j = pd.Series(j, name="j", dtype=str)
        values = pd.Series(np.hstack(values), name="value")

        return pd.concat([i, j, values], axis=1)

    def plot(self, cmap=None, title=""):
        """Create a heatmap of the dissimilarity matrix.

        Parameters
        ----------
        cmap: str or matplotlib.colors.Colormap, optional
            Sets the color scheme of the heatmap
            If ``None``, defaults to the colormap specified in the matplotlib
            rc file.

        title: str, optional
            Sets the title label of the heatmap
            (Default is blank)

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing the heatmap and colorbar of the plotted
            dissimilarity matrix.

        Examples
        --------
        .. plot::

           Define a dissimilarity matrix with five objects labeled A-E:

           >>> from skbio.stats.distance import DissimilarityMatrix
           >>> dm = DissimilarityMatrix([[0, 1, 2, 3, 4], [1, 0, 1, 2, 3],
           ...                           [2, 1, 0, 1, 2], [3, 2, 1, 0, 1],
           ...                           [4, 3, 2, 1, 0]],
           ...                          ['A', 'B', 'C', 'D', 'E'])

           Plot the dissimilarity matrix as a heatmap:

           >>> fig = dm.plot(cmap='Reds', title='Example heatmap')  # doctest: +SKIP

        """
        self._get_mpl_plt()

        # based on http://stackoverflow.com/q/14391959/3776794
        fig, ax = self.plt.subplots()

        # use pcolormesh instead of pcolor for performance
        heatmap = ax.pcolormesh(self.data, cmap=cmap)
        fig.colorbar(heatmap)

        # center labels within each cell
        ticks = np.arange(0.5, self.shape[0])
        ax.set_xticks(ticks, minor=False)
        ax.set_yticks(ticks, minor=False)

        # Ensure there is no white border around the heatmap by manually
        # setting the limits
        ax.set_ylim(0, len(self.ids))
        ax.set_xlim(0, len(self.ids))

        # display data as it is stored in the dissimilarity matrix
        # (default is to have y-axis inverted)
        ax.invert_yaxis()

        ax.set_xticklabels(self.ids, rotation=90, minor=False)
        ax.set_yticklabels(self.ids, minor=False)

        ax.set_title(title)

        return fig

    def to_data_frame(self):
        """Create a ``pandas.DataFrame`` from this ``DissimilarityMatrix``.

        Returns
        -------
        pd.DataFrame
            ``pd.DataFrame`` with IDs on index and columns.

        Examples
        --------
        >>> from skbio import DistanceMatrix
        >>> dm = DistanceMatrix([[0, 1, 2],
        ...                      [1, 0, 3],
        ...                      [2, 3, 0]], ids=['a', 'b', 'c'])
        >>> df = dm.to_data_frame()
        >>> df
             a    b    c
        a  0.0  1.0  2.0
        b  1.0  0.0  3.0
        c  2.0  3.0  0.0

        """
        return pd.DataFrame(data=self.data, index=self.ids, columns=self.ids)

    def __str__(self):
        """Return a string representation of the dissimilarity matrix.

        Summary includes matrix dimensions, a (truncated) list of IDs, and
        (truncated) array of dissimilarities.

        Returns
        -------
        str
            String representation of the dissimilarity matrix.

        """
        return "%dx%d %s matrix\nIDs:\n%s\nData:\n" % (
            self.shape[0],
            self.shape[1],
            self._matrix_element_name,
            _pprint_strs(self.ids),
        ) + str(self.data)

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

        """
        if isinstance(index, str):
            return self.data[self.index(index)]
        elif self._is_id_pair(index):
            return self.data[self.index(index[0]), self.index(index[1])]
        else:
            return self.data.__getitem__(index)

    def _validate_ids(self, data, ids):
        """Validate the IDs.

        Checks that IDs are unique and that the number of IDs matches the
        number of rows/cols in the data array.

        Subclasses can override this method to perform different/more specific
        validation.

        Notes
        -----
        Accepts arguments instead of inspecting instance attributes to avoid
        creating an invalid dissimilarity matrix before raising an error.
        Otherwise, the invalid dissimilarity matrix could be used after the
        exception is caught and handled.

        """
        duplicates = find_duplicates(ids)
        if duplicates:
            formatted_duplicates = ", ".join(repr(e) for e in duplicates)
            raise DissimilarityMatrixError(
                "IDs must be unique. Found the "
                "following duplicate IDs: %s" % formatted_duplicates
            )
        if 0 == len(ids):
            raise DissimilarityMatrixError("IDs must be at least 1 in " "size.")
        if len(ids) != data.shape[0]:
            raise DissimilarityMatrixError(
                "The number of IDs (%d) must match "
                "the number of rows/columns in the "
                "data (%d)." % (len(ids), data.shape[0])
            )

    def _validate_shape(self, data):
        """Validate the data array shape.

        Checks that the data is at least 1x1 in size, 2D, square, and
        contains only floats.

        Notes
        -----
        Accepts arguments instead of inspecting instance attributes to avoid
        creating an invalid dissimilarity matrix before raising an error.
        Otherwise, the invalid dissimilarity matrix could be used after the
        exception is caught and handled.

        """
        if 0 in data.shape:
            raise DissimilarityMatrixError("Data must be at least 1x1 in " "size.")
        if len(data.shape) != 2:
            raise DissimilarityMatrixError("Data must have exactly two " "dimensions.")
        if data.shape[0] != data.shape[1]:
            raise DissimilarityMatrixError(
                "Data must be square (i.e., have "
                "the same number of rows and "
                "columns)."
            )
        if data.dtype not in (np.float32, np.float64):
            raise DissimilarityMatrixError(
                "Data must contain only floating " "point values."
            )

    def _validate(self, data, ids):
        """Validate the data array and IDs.

        Checks that the data is at least 1x1 in size, 2D, square, and
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
        self._validate_shape(data)
        self._validate_ids(data, ids)

    def _index_list(self, list_):
        return {id_: idx for idx, id_ in enumerate(list_)}

    def _is_id_pair(self, index):
        return (
            isinstance(index, tuple)
            and len(index) == 2
            and all(map(lambda e: isinstance(e, str), index))
        )


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
    _matrix_element_name = "distance"

    @classonlymethod
    def from_iterable(cls, iterable, metric, key=None, keys=None, validate=True):
        """Create DistanceMatrix from all pairs in an iterable given a metric.

        Parameters
        ----------
        iterable : iterable
            Iterable containing objects to compute pairwise distances on.
        metric : callable
            A function that takes two arguments and returns a float
            representing the distance between the two arguments.
        key : callable or metadata key, optional
            A function that takes one argument and returns a string
            representing the id of the element in the distance matrix.
            Alternatively, a key to a `metadata` property if it exists for
            each element in the `iterable`. If None, then default ids will be
            used.
        keys : iterable, optional
            An iterable of the same length as `iterable`. Each element will be
            used as the respective key.
        validate : boolean, optional
            If ``True``, all pairwise distances are computed, including upper
            and lower triangles and the diagonal, and the resulting matrix is
            validated for symmetry and hollowness. If ``False``, `metric` is
            assumed to be hollow and symmetric and only the lower triangle
            (excluding the diagonal) is computed. Pass ``validate=False`` if
            you are sure `metric` is hollow and symmetric for improved
            performance.

        Returns
        -------
        DistanceMatrix
            The `metric` applied to pairwise elements in the `iterable`.

        Raises
        ------
        ValueError
            If `key` and `keys` are both provided.

        """
        if validate:
            return super(DistanceMatrix, cls).from_iterable(iterable, metric, key, keys)

        iterable = list(iterable)
        if key is not None and keys is not None:
            raise ValueError("Cannot use both `key` and `keys` at the same" " time.")

        keys_ = None
        if key is not None:
            keys_ = [resolve_key(e, key) for e in iterable]
        elif keys is not None:
            keys_ = keys

        dm = np.zeros((len(iterable),) * 2)
        for i, a in enumerate(iterable):
            for j, b in enumerate(iterable[:i]):
                dm[i, j] = dm[j, i] = metric(a, b)

        return cls(dm, keys_)

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
        return squareform(self._data, force="tovector", checks=False)

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

        if condensed:
            permuted_condensed = distmat_reorder_condensed(self._data, order)
            return permuted_condensed
        else:
            # Note: Skip validation, since we assume self was already validated
            permuted = distmat_reorder(self._data, order)
            return self.__class__(permuted, self.ids, validate=False)

    def _validate(self, data, ids):
        """Validate the data array and IDs.

        Overrides the superclass `_validate`. Performs a check for symmetry in
        addition to the checks performed in the superclass.

        """
        super(DistanceMatrix, self)._validate(data, ids)

        data_sym, data_hol = is_symmetric_and_hollow(data)

        if not data_sym:
            raise DistanceMatrixError("Data must be symmetric and cannot contain NaNs.")

        if not data_hol:
            raise DistanceMatrixError(
                "Data must be hollow (i.e., the diagonal" " can only contain zeros)."
            )

    def to_series(self):
        """Create a ``pandas.Series`` from this ``DistanceMatrix``.

        The series will contain distances in condensed form: only distances
        from one matrix triangle are included, and the diagonal is excluded.
        The series' index will be a ``pd.MultiIndex`` relating pairs of IDs to
        distances. The pairs of IDs will be in row-major order with respect to
        the upper matrix triangle.

        To obtain all distances (i.e. both upper and lower matrix triangles and
        the diagonal), use ``DistanceMatrix.to_data_frame``. To obtain *only*
        the distances in condensed form (e.g. for use with SciPy), use
        ``DistanceMatrix.condensed_form``.

        Returns
        -------
        pd.Series
            ``pd.Series`` with pairs of IDs on the index.

        See Also
        --------
        to_data_frame
        condensed_form
        scipy.spatial.distance.squareform

        Examples
        --------
        >>> from skbio import DistanceMatrix
        >>> dm = DistanceMatrix([[0, 1, 2, 3],
        ...                      [1, 0, 4, 5],
        ...                      [2, 4, 0, 6],
        ...                      [3, 5, 6, 0]], ids=['a', 'b', 'c', 'd'])
        >>> dm.to_series()
        a  b    1.0
           c    2.0
           d    3.0
        b  c    4.0
           d    5.0
        c  d    6.0
        dtype: float64

        """
        distances = self.condensed_form()
        # `id_pairs` will not be interpreted as a `pd.MultiIndex` if it is an
        # iterable returned by `itertools.combinations`.
        id_pairs = list(itertools.combinations(self.ids, 2))
        index = pd.Index(id_pairs, tupleize_cols=True)
        return pd.Series(data=distances, index=index, dtype=float)


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


# helper functions for anosim and permanova


def _preprocess_input_sng(ids, sample_size, grouping, column):
    """Compute intermediate results not affected by permutations.

    These intermediate results can be computed a single time for efficiency,
    regardless of grouping vector permutations (i.e., when calculating the
    p-value). These intermediate results are used by both ANOSIM and PERMANOVA.

    Also validates and normalizes input (e.g., converting ``DataFrame`` column
    into grouping vector).

    """
    if isinstance(grouping, pd.DataFrame):
        if column is None:
            raise ValueError("Must provide a column name if supplying a DataFrame.")
        else:
            grouping = _df_to_vector(ids, grouping, column)
    elif isinstance(grouping, pd.Series):
        if (column is not None) and (column != grouping.name):
            raise ValueError(
                "Column name does not match your Series name. Try not"
                " providing column at all."
            )
        else:
            grouping = _df_to_vector(ids, grouping.to_frame(), column=grouping.name)
    elif column is not None:
        raise ValueError("Must provide a DataFrame if supplying a column name.")

    if len(grouping) != sample_size:
        raise ValueError(
            "Grouping vector size must match the number of IDs in the "
            "distance matrix."
        )

    # Find the group labels and convert grouping to an integer vector
    # (factor).
    groups, grouping = np.unique(grouping, return_inverse=True)
    num_groups = len(groups)

    if num_groups == len(grouping):
        raise ValueError(
            "All values in the grouping vector are unique. This method cannot "
            "operate on a grouping vector with only unique values (e.g., "
            "there are no 'within' distances because each group of objects "
            "contains only a single object)."
        )
    if num_groups == 1:
        raise ValueError(
            "All values in the grouping vector are the same. This method "
            "cannot operate on a grouping vector with only a single group of "
            "objects (e.g., there are no 'between' distances because there is "
            "only a single group)."
        )

    return num_groups, grouping


def _preprocess_input(distance_matrix, grouping, column):
    """Compute intermediate results not affected by permutations.

    These intermediate results can be computed a single time for efficiency,
    regardless of grouping vector permutations (i.e., when calculating the
    p-value). These intermediate results are used by both ANOSIM and PERMANOVA.

    Also validates and normalizes input (e.g., converting ``DataFrame`` column
    into grouping vector).

    """
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("Input must be a DistanceMatrix.")
    sample_size = distance_matrix.shape[0]

    num_groups, grouping = _preprocess_input_sng(
        distance_matrix.ids, sample_size, grouping, column
    )

    tri_idxs = np.triu_indices(sample_size, k=1)
    distances = distance_matrix.condensed_form()

    return sample_size, num_groups, grouping, tri_idxs, distances


def _df_to_vector(ids, df, column):
    """Return a grouping vector from a ``DataFrame`` column.

    Parameters
    ----------
    ids : liat
        IDs that will be mapped to group labels.
    df : pandas.DataFrame
        ``DataFrame`` (indexed by distance matrix ID).
    column : str
        Column name in `df` containing group labels.

    Returns
    -------
    list
        Grouping vector (vector of labels) based on the IDs in
        `ids`. Each ID's label is looked up in the ``DataFrame``
        under the column specified by `column`.

    Raises
    ------
    ValueError
        If `column` is not in the ``DataFrame``, or a distance matrix ID is
        not in the ``DataFrame``.

    """
    if column not in df:
        raise ValueError("Column '%s' not in DataFrame." % column)

    grouping = df.reindex(ids, axis=0).loc[:, column]
    if grouping.isnull().any():
        raise ValueError(
            "One or more IDs in the distance matrix are not in the data " "frame."
        )
    return grouping.tolist()


def _run_monte_carlo_stats(test_stat_function, grouping, permutations):
    """Run stat test and compute significance with Monte Carlo permutations."""
    if permutations < 0:
        raise ValueError(
            "Number of permutations must be greater than or equal to zero."
        )

    stat = test_stat_function(grouping)

    p_value = np.nan
    if permutations > 0:
        perm_stats = np.empty(permutations, dtype=np.float64)

        for i in range(permutations):
            perm_grouping = np.random.permutation(grouping)
            perm_stats[i] = test_stat_function(perm_grouping)

        p_value = ((perm_stats >= stat).sum() + 1) / (permutations + 1)

    return stat, p_value


def _build_results(
    method_name, test_stat_name, sample_size, num_groups, stat, p_value, permutations
):
    """Return ``pandas.Series`` containing results of statistical test."""
    return pd.Series(
        data=[
            method_name,
            test_stat_name,
            sample_size,
            num_groups,
            stat,
            p_value,
            permutations,
        ],
        index=[
            "method name",
            "test statistic name",
            "sample size",
            "number of groups",
            "test statistic",
            "p-value",
            "number of permutations",
        ],
        name="%s results" % method_name,
    )
