# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
from copy import deepcopy
from typing import (
    Any,
    Callable,
    Iterable,
    ClassVar,
    Collection,
    Type,
    Optional,
    Sequence,
    Union,
    TYPE_CHECKING,
)

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform

from skbio._base import SkbioObject
from skbio.stats._misc import _pprint_strs
from skbio.util import find_duplicates, get_rng
from skbio.util._decorator import classonlymethod
from skbio.util._misc import resolve_key
from skbio.util._plotting import PlottableMixin
from skbio.io.descriptors import Read, Write

from ._utils import is_symmetric_and_hollow, is_symmetric
from ._utils import distmat_reorder, distmat_reorder_condensed

if TYPE_CHECKING:  # pragma: no cover
    from numpy.random import Generator
    import matplotlib.figure
    from matplotlib.colors import Colormap
    from skbio.util._typing import SeedLike


class PairwiseMatrixError(Exception):
    """General error for dissimilarity matrix validation failures."""

    pass


class SymmetricMatrixError(PairwiseMatrixError):
    """General error for symmetric matrix validation failures."""

    pass


class DistanceMatrixError(SymmetricMatrixError):
    """General error for distance matrix validation failures."""

    pass


class MissingIDError(PairwiseMatrixError):
    """Error for ID lookup that doesn't exist in the pairwise matrix."""

    def __init__(self, missing_id):
        super(MissingIDError, self).__init__()
        self.args = ("The ID '%s' is not in the matrix." % missing_id,)


class PairwiseMatrix(SkbioObject, PlottableMixin):
    r"""Store pairwise relationships between objects.

    A `PairwiseMatrix` instance stores a square, two-dimensional
    matrix of relationships between objects. Objects could be, for example,
    samples or DNA sequences. A sequence of IDs accompanies the
    data.

    Methods are provided to load and save pairwise matrices from/to disk,
    as well as perform common operations such as extracting values
    based on object ID. Additionally, the
    :meth:`plot` method provides
    convenient built-in plotting functionality.

    Parameters
    ----------
    data : array_like or PairwiseMatrix
        Square, two-dimensional ``numpy.ndarray`` of pairwise relationships
        (floats), or a structure that can be converted to a ``numpy.ndarray``
        using ``numpy.asarray`` or a one-dimensional vector of pairwise relationships
        (floats), as defined by `scipy.spatial.distance.squareform
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html>`_.
        Can instead be a `PairwiseMatrix` (or subclass) instance, in which case the
        instance's data will be used. Data will be converted to a float ``dtype`` if
        necessary. A copy will *not* be made if already a ``numpy.ndarray`` with a
        float ``dtype``.
    ids : sequence of str, optional
        Sequence of strings to be used as object IDs. Must match the number of
        rows/cols in `data`. If ``None`` (default), IDs will be
        monotonically-increasing integers cast as strings, with numbering
        starting from zero, e.g., ``('0', '1', '2', '3', ...)``.
    validate : bool, optional
        If `validate` is ``True`` (default) and data is not a
        `PairwiseMatrix` object, the input data will be validated.

    See Also
    --------
    SymmetricMatrix
    DistanceMatrix
    scipy.spatial.distance.squareform

    Notes
    -----
    The values are stored in redundant (square-form) format [1]_.

    The data are not checked for symmetry or hollowness, nor guaranteed/assumed to be
    symmetric or hollow.

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    """

    default_write_format: ClassVar[str] = "lsmat"
    # Used in __str__
    # TODO: decide on what to call a matrix element here
    _matrix_element_name: ClassVar[str] = "dissimilarity"

    read = Read()
    write = Write()

    def __init__(
        self,
        data: Union[
            np.ndarray,
            Sequence[float],
            Sequence[Sequence[float]],
            "PairwiseMatrix",
        ],
        ids: Optional[Sequence[str]] = None,
        validate: bool = True,
    ) -> None:
        data, ids = self._normalize_input(data, ids)
        # convert data to redundant if 1D input.
        # should do this for PairwiseMatrix only.
        if data.ndim == 1:
            data = squareform(data, force="tomatrix", checks=False)

        if ids is None:
            ids = self._generate_ids(data)
        else:
            ids = tuple(ids)

        # validation
        if validate:
            self._validate_shape(data)
            self._validate_ids(data, ids)

        self._ids = ids
        self._id_index = self._index_list(self._ids)
        self._data = self._init_data(data)
        self._flags = self._init_flags()

    def _normalize_input(self, data, ids):
        """Get input into standard numpy array format."""
        if isinstance(data, PairwiseMatrix):
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

        # Make data_ explicitly an ndarray to help with type checking.
        # At this point in the code we can be certain that data is an
        # ndarray.
        assert isinstance(data, np.ndarray)
        data_: np.ndarray = data
        return data_, ids

    def _generate_ids(self, data):
        """Generate ids if none provided."""
        if data.ndim != 1:
            return tuple(str(i) for i in range(data.shape[0]))
        else:
            return tuple(str(i) for i in range(_vec_to_shape(data)))

    def _init_flags(self):
        """Initialize boolean flags for matrix forms."""
        # This is the default value. PairwiseMatrix doesn't really need flags.
        return {"CONDENSED": False}

    def _init_data(self, data):
        """Initialize underlying data structure."""
        return data

    @classonlymethod
    def from_iterable(
        cls,
        iterable: Iterable[Any],
        metric: Callable,
        key: Optional[Any] = None,
        keys: Optional[Iterable[Any]] = None,
    ) -> "PairwiseMatrix":
        r"""Create PairwiseMatrix from an iterable given a metric.

        Parameters
        ----------
        iterable : iterable
            Iterable containing objects to compute pairwise relationships on.
        metric : callable
            A function that takes two arguments and returns a float
            representing the relationship between the two arguments.
        key : callable or metadata key, optional
            A function that takes one argument and returns a string
            representing the id of the element in the pairwise matrix.
            Alternatively, a key to a `metadata` property if it exists for
            each element in the `iterable`. If None, then default ids will be
            used.
        keys : iterable, optional
            An iterable of the same length as `iterable`. Each element will be
            used as the respective key.

        Returns
        -------
        PairwiseMatrix
            The `metric` applied to all pairwise elements in the `iterable`.

        Raises
        ------
        ValueError
            If `key` and `keys` are both provided.

        """
        iterable = list(iterable)
        if key is not None and keys is not None:
            raise ValueError("Cannot use both `key` and `keys` at the same time.")

        keys_ = None
        if key is not None:
            keys_ = [resolve_key(e, key) for e in iterable]
        elif keys is not None:
            keys_ = list(keys)

        dm = np.empty((len(iterable),) * 2)
        for i, a in enumerate(iterable):
            for j, b in enumerate(iterable):
                dm[i, j] = metric(a, b)

        return cls(dm, keys_)  # type: ignore[operator]

    @property
    def data(self) -> np.ndarray:
        r"""Array of pairwise relationships.

        A square, two-dimensional ``numpy.ndarray`` of values (floats). A copy is *not*
        returned.

        Notes
        -----
        This property is not writeable.

        """
        return self._data

    @property
    def ids(self) -> tuple:
        r"""Tuple of object IDs.

        A tuple of strings, one for each object in the pairwise matrix.

        Notes
        -----
        This property is writeable, but the number of new IDs must match the
        number of objects in `data`.

        """
        return self._ids

    @ids.setter
    def ids(self, ids_: Sequence[str]) -> None:
        ids_ = tuple(ids_)
        self._validate_ids(self._data, ids_)
        self._ids = ids_
        self._id_index = self._index_list(self._ids)

    @property
    def dtype(self) -> np.dtype:
        """Data type of the matrix values."""
        return self._data.dtype

    @property
    def shape(self) -> tuple:
        r"""Two-element tuple containing the redundant form matrix dimensions.

        Notes
        -----
        As the matrix is guaranteed to be square, both tuple entries will always be
        equal. The shape of the redundant form matrix is returned.

        """
        if self._flags["CONDENSED"]:
            m = _vec_to_shape(self._data)
            return (m, m)
        return self._data.shape

    @property
    def size(self) -> int:
        r"""Total number of elements in the underlying data structure.

        Notes
        -----
        If the matrix is stored in redundant form, size is equivalent to
        ``self.shape[0] * self.shape[1]``. If the matrix is stored in condensed form,
        size is equal to the number of elements in the condensed array.

        """
        return self._data.size

    @property
    def T(self) -> "PairwiseMatrix":
        r"""Transpose of the matrix.

        See Also
        --------
        transpose

        """
        return self.transpose()

    def transpose(self) -> "PairwiseMatrix":
        r"""Return the transpose of the matrix.

        Notes
        -----
        A deep copy is returned.

        Returns
        -------
        PairwiseMatrix
            Transpose of the matrix. Will be the same type as ``self``.

        """
        # Note: Skip validation, since we assume self was already validated
        return self._copy(transpose=True)

    def index(self, lookup_id: str) -> int:
        r"""Return the index of the specified ID.

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
            If `lookup_id` is not in the matrix.

        """
        if lookup_id in self:
            return self._id_index[lookup_id]
        else:
            raise MissingIDError(lookup_id)

    def redundant_form(self) -> np.ndarray:
        r"""Return an array of values in redundant form.

        Returns
        -------
        ndarray
            Two-dimensional `numpy.ndarray` of values in redundant
            form.

        Notes
        -----
        Redundant form is described in [1]_.

        Does *not* return a copy of the data.

        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

        """
        return self._data

    def copy(self) -> "PairwiseMatrix":
        r"""Return a deep copy of the matrix.

        Returns
        -------
        PairwiseMatrix
            Deep copy of the matrix. Will be the same type as
            `self`.

        """
        # We deepcopy IDs in case the tuple contains mutable objects at some
        # point in the future.
        # Note: Skip validation, since we assume self was already validated
        # We deepcopy IDs in case the tuple contains mutable objects at some
        # point in the future.
        # Note: Skip validation, since we assume self was already validated
        return self._copy()

    def _copy(self, transpose: bool = False) -> "PairwiseMatrix":
        """Copy support.

        Parameters
        ----------
        transpose : bool
            Transpose the data on copy.

        Returns
        -------
        PairwiseMatrix
            Deep copy of the matrix. Will be the same type as `self`.
        """
        data = self._data.copy()
        if transpose:
            data = data.T
        return self.__class__(data, deepcopy(self.ids), validate=False)

    def rename(self, mapper: Union[dict, Callable], strict: bool = True) -> None:
        r"""Rename IDs in the matrix.

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
        >>> from skbio.stats.distance import PairwiseMatrix
        >>> dm = PairwiseMatrix([[0, 1], [2, 3]], ids=['a', 'b'])
        >>> dm.rename({'a': 'x', 'b': 'y'})
        >>> print(dm.ids)
        ('x', 'y')

        """
        if isinstance(mapper, dict):
            if strict and not set(self.ids).issubset(mapper):
                raise ValueError(
                    "The IDs in mapper do not include all IDs in the matrix."
                )
            new_ids = tuple(mapper.get(x, x) for x in self.ids)
        else:
            new_ids = tuple(mapper(x) for x in self.ids)
        self.ids = new_ids

    def filter(
        self,
        ids: Sequence[str],
        strict: bool = True,  # , preserve_condensed=True
    ) -> "PairwiseMatrix":
        r"""Filter the matrix by IDs.

        Parameters
        ----------
        ids : sequence of str
            IDs to retain. May not contain duplicates or be empty. Each ID must
            be present in the matrix.
        strict : bool, optional
            If `strict` is ``True`` and an ID that is not found in the distance
            matrix is found in `ids`, a ``MissingIDError`` exception will be
            raised, otherwise the ID will be ignored.

        Returns
        -------
        PairwiseMatrix
            Filtered matrix containing only the IDs specified in
            `ids`. IDs will be in the same order as they appear in `ids`.

        Raises
        ------
        MissingIDError
            If an ID in `ids` is not in the object's list of IDs.

        """
        if tuple(self._ids) == tuple(ids):
            if self._flags["CONDENSED"]:
                return self.__class__(self._data, self._ids, condensed=True)
            else:
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
        if self._flags["CONDENSED"]:
            filtered_data = distmat_reorder_condensed_py(self.condensed_form(), idxs)
            self._validate_ids(filtered_data, ids)
            return self.__class__(filtered_data, ids, validate=False, condensed=True)
        else:
            filtered_data = distmat_reorder(self.redundant_form(), idxs)
            self._validate_ids(filtered_data, ids)
            return self.__class__(filtered_data, ids, validate=False)

    def _stable_order(self, ids: Iterable[str]) -> np.ndarray:
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

    def within(self, ids: Iterable[str]) -> pd.DataFrame:
        r"""Obtain all the pairwise values among the set of IDs.

        Parameters
        ----------
        ids : iterable of str
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
            If an ID(s) specified is not in the matrix.

        Notes
        -----
        Order of the return items is stable, meaning that requesting IDs
        ['a', 'b'] is equivalent to ['b', 'a']. The order is with respect
        to the order of the .ids attribute of self.

        Examples
        --------
        >>> from skbio.stats.distance import PairwiseMatrix
        >>> dm = PairwiseMatrix([[0, 1, 2, 3, 4], [1, 0, 1, 2, 3],
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
                "At least one ID (e.g., '%s') was not found." % not_present.pop()
            )

        return self._subset_to_dataframe(ids, ids)

    def between(
        self, from_: Iterable[str], to_: Iterable[str], allow_overlap: bool = False
    ) -> pd.DataFrame:
        """Obtain the pairwise values between the two groups of IDs.

        Parameters
        ----------
        from_ : Iterable of str
            The IDs to obtain values from. Values from all pairs of IDs
            in from and to will be obtained.
        to_ : Iterable of str
            The IDs to obtain values to. Values from all pairs of IDs
            in to and from will be obtained.
        allow_overlap : bool, optional
            If True, allow overlap in the IDs of from and to (which would in
            effect be collecting the within distances). Default is False.

        Returns
        -------
        pd.DataFrame
            (i, j, value) representing the source ID ("i"), the target ID ("j")
            and the value between them ("value").

        Raises
        ------
        MissingIDError
            If an ID(s) specified is not in the matrix.

        Notes
        -----
        Order of the return items is stable, meaning that requesting IDs
        ['a', 'b'] is equivalent to ['b', 'a']. The order is with respect to
        the .ids attribute of self.

        Examples
        --------
        >>> from skbio.stats.distance import PairwiseMatrix
        >>> dm = PairwiseMatrix([[0, 1, 2, 3, 4], [1, 0, 1, 2, 3],
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
                "At least one ID (e.g., '%s') was not found." % not_present.pop()
            )

        overlapping = from_ & to_
        if not allow_overlap and overlapping:
            raise KeyError(
                "At least one ID overlaps in from_ and to_ "
                "(e.g., '%s'). This constraint can removed with "
                "allow_overlap=True." % overlapping.pop()
            )

        return self._subset_to_dataframe(from_, to_)

    def _subset_to_dataframe(
        self, i_ids: Iterable[str], j_ids: Iterable[str]
    ) -> pd.DataFrame:
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
        j_labels = tuple(self.ids[j] for j in j_indices)

        i: list[str] = []
        j: list[str] = []

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

    def plot(
        self, cmap: Optional[Union[str, "Colormap"]] = None, title: str = ""
    ) -> "matplotlib.figure.Figure":
        """Create a heatmap of the matrix.

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
            Figure containing the heatmap and colorbar of the plotted matrix.

        Examples
        --------
        .. plot::

           Define a matrix with five objects labeled A-E:

           >>> from skbio.stats.distance import PairwiseMatrix
           >>> dm = PairwiseMatrix([[0, 1, 2, 3, 4], [1, 0, 1, 2, 3],
           ...                           [2, 1, 0, 1, 2], [3, 2, 1, 0, 1],
           ...                           [4, 3, 2, 1, 0]],
           ...                          ['A', 'B', 'C', 'D', 'E'])

           Plot the matrix as a heatmap:

           >>> fig = dm.plot(cmap='Reds', title='Example heatmap')  # doctest: +SKIP

        """
        self._get_mpl_plt()

        # based on http://stackoverflow.com/q/14391959/3776794
        fig, ax = self.plt.subplots()

        # use pcolormesh instead of pcolor for performance
        heatmap = ax.pcolormesh(self.redundant_form(), cmap=cmap)
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

    def to_data_frame(self) -> pd.DataFrame:
        """Create a ``pandas.DataFrame`` from this ``PairwiseMatrix``.

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
        return pd.DataFrame(
            data=self.redundant_form(), index=self.ids, columns=self.ids
        )

    def __str__(self) -> str:
        """Return a string representation of the matrix.

        Summary includes matrix dimensions, a (truncated) list of IDs, and
        (truncated) array of values.

        Returns
        -------
        str
            String representation of the matrix.

        """
        return "%dx%d %s matrix\nIDs:\n%s\nData:\n" % (
            self.shape[0],
            self.shape[1],
            self._matrix_element_name,
            _pprint_strs(self.ids),
        ) + str(self._data)

    def __eq__(self, other: object) -> bool:
        r"""Compare this matrix to another for equality.

        Two matrices are equal if they have the same shape, IDs
        (in the same order!), and have data arrays that are equal.

        Checks are *not* performed to ensure that `other` is a
        `PairwiseMatrix` instance.

        Parameters
        ----------
        other : PairwiseMatrix
            Matrix to compare to for equality.

        Returns
        -------
        bool
            ``True`` if `self` is equal to `other`, ``False`` otherwise.

        """
        if not isinstance(other, PairwiseMatrix):
            return NotImplemented

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
            elif not np.array_equal(self._data, other.data):
                equal = False
        except AttributeError:
            equal = False

        return equal

    def __ne__(self, other: object) -> bool:
        """Determine whether two matrices are not equal.

        Parameters
        ----------
        other : PairwiseMatrix
            Matrix to compare to.

        Returns
        -------
        bool
            ``True`` if `self` is not equal to `other`, ``False`` otherwise.

        See Also
        --------
        __eq__

        """
        return not self == other

    def __contains__(self, lookup_id: str) -> bool:
        r"""Check if the specified ID is in the matrix.

        Parameters
        ----------
        lookup_id : str
            ID to search for.

        Returns
        -------
        bool
            `True` if `lookup_id` is in the matrix, `False` otherwise.

        See Also
        --------
        index

        """
        return lookup_id in self._id_index

    def __getitem__(
        self, index: Union[str, tuple[str, str], Any]
    ) -> Union[np.ndarray, float]:
        r"""Slice into data by object ID or numpy indexing.

        Extracts data from the matrix by object ID, a pair of
        IDs, or numpy indexing/slicing.

        Parameters
        ----------
        index : str, two-tuple of str, or numpy index
            `index` can be one of the following forms: an ID, a pair of IDs, or
            a numpy index.

            If `index` is a string, it is assumed to be an ID and a
            ``numpy.ndarray`` row vector is returned for the corresponding ID.
            Note that the ID's row of values is returned, *not* its
            column. If the matrix is symmetric, the two will be identical, but
            this makes a difference if the matrix is asymmetric.

            If `index` is a two-tuple of strings, each string is assumed to be
            an ID and the corresponding matrix element is returned that
            represents the value between the two IDs. Note that the
            order of lookup by ID pair matters if the matrix is asymmetric: the
            first ID will be used to look up the row, and the second ID will be
            used to look up the column. Thus, ``dm['a', 'b']`` may not be the
            same as ``dm['b', 'a']`` if the matrix is asymmetric.

            Otherwise, `index` will be passed through to
            ``PairwiseMatrix.data.__getitem__``, allowing for standard
            indexing of a ``numpy.ndarray`` (e.g., slicing).

        Returns
        -------
        ndarray or scalar
            Indexed data, where return type depends on the form of `index` (see
            description of `index` for more details).

        Raises
        ------
        MissingIDError
            If the ID(s) specified in `index` are not in the matrix.

        Notes
        -----
        The lookup based on ID(s) is quick. NumPy indexing (slicing) on condensed
        form matrices will convert them to redundant, roughly doubling their memory
        requirement.

        """
        if isinstance(index, str):
            row_idx = self.index(index)

            return self._data[row_idx]
        elif isinstance(index, tuple) and self._is_id_pair(index):
            i, j = self.index(index[0]), self.index(index[1])
            return self._data[i, j]
        else:
            # NumPy index types are numerous and complex, easier to just
            # ignore them in type checking.
            return self._data.__getitem__(index)  # type: ignore[index]

    def _validate_ids(self, data: np.ndarray, ids: Collection[str]) -> None:
        """Validate the IDs.

        Checks that IDs are unique and that the number of IDs matches the
        number of rows/cols in the data array.

        Subclasses can override this method to perform different/more specific
        validation.

        Notes
        -----
        Accepts arguments instead of inspecting instance attributes to avoid
        creating an invalid matrix before raising an error. Otherwise, the invalid
        matrix could be used after the exception is caught and handled.

        """
        duplicates = find_duplicates(ids)
        if duplicates:
            formatted_duplicates = ", ".join(repr(e) for e in duplicates)
            raise PairwiseMatrixError(
                "IDs must be unique. Found the "
                "following duplicate IDs: %s" % formatted_duplicates
            )
        if 0 == len(ids):
            raise PairwiseMatrixError("IDs must be at least 1 in size.")
        # handle condensed form data
        if data.ndim == 1:
            n = _vec_to_shape(data)
            if len(ids) != n:
                raise PairwiseMatrixError(
                    "The number of IDs (%d) must match "
                    "the number of rows/columns in the "
                    "data (%d)." % (len(ids), n)
                )
        # handle redundant form data
        else:
            if len(ids) != data.shape[0]:
                raise PairwiseMatrixError(
                    "The number of IDs (%d) must match "
                    "the number of rows/columns in the "
                    "data (%d)." % (len(ids), data.shape[0])
                )

    def _validate_shape(self, data: np.ndarray) -> None:
        """Validate the data array shape.

        Checks that the data is at least 1x1 in size, 2D, square, and
        contains only floats.

        Notes
        -----
        Accepts arguments instead of inspecting instance attributes to avoid
        creating an invalid matrix before raising an error. Otherwise, the invalid
        matrix could be used after the exception is caught and handled.

        """
        if 0 in data.shape:
            raise PairwiseMatrixError("Data must be at least 1x1 in size.")
        if len(data.shape) != 2:
            raise PairwiseMatrixError("Data must have exactly two dimensions.")
        if data.shape[0] != data.shape[1]:
            raise PairwiseMatrixError(
                "Data must be square (i.e., have the same number of rows and columns)."
            )
        if data.dtype not in (np.float32, np.float64):
            raise PairwiseMatrixError("Data must contain only floating point values.")

    def _index_list(self, list_: Sequence[str]) -> dict:
        return {id_: idx for idx, id_ in enumerate(list_)}

    def _is_id_pair(self, index: tuple) -> bool:
        return len(index) == 2 and all(map(lambda e: isinstance(e, str), index))


class SymmetricMatrix(PairwiseMatrix):
    """Store symmetric pairwise relationships between objects.

    A `SymmetricMatrix` is a `PairwiseMatrix` with the additional requirement that the
    matrix data is symmetric. There are additional methods made available that take
    advantage of this symmetry. The :func:`~skbio.stats.distance.PairwiseMatrix.plot`
    method provides convenient built-in plotting functionality.

    Parameters
    ----------
    data : array_like or SymmetricMatrix
        Square, symmetric, two-dimensional ``numpy.ndarray`` values (floats), or a
        structure that can be converted to a ``numpy.ndarray`` using ``numpy.asarray``
        or a one-dimensional vector of values (floats), as defined by
        `scipy.spatial.distance.squareform
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squaref>`_.
        Can instead be a `PairwiseMatrix` (or `SymmetricMatrix`) instance, in which the
        instance's data will be used. Data will be converted to a float ``dtype`` if
        necessary. A copy will *not* be made if already a ``numpy.ndarray`` with a
        float ``dtype``.
    ids : sequence of str, optional
        Sequence of strings to be used as object IDs. Must match the number of
        rows/cols in `data`. If ``None`` (the default), IDs will be
        monotonically-increasing integers cast as strings, with numbering
        starting from zero, e.g., ``('0', '1', '2', '3', ...)``.
    validate : bool, optional
        If `validate` is ``True`` (the default) and data is not a
        DistanceMatrix object, the input data will be validated.
    condensed : bool, optional
        If ``False`` (default) the data will be stored in a 2-dimensional redundant
        form. If ``True``, the data will be stored in a 1-dimensional condensed form.
    diagonal : np.ndarray or float, optional
        If the matrix is to be stored in condensed form, `diagonal` stores the
        information contained in the matrix's diagonal. If `diagonal` is a float, this
        value represents all values in the diagonal of the matrix. If `diagonal` is an
        array, it's values represent the values along the diagonal of the original
        matrix.

    See Also
    --------
    PairwiseMatrix
    SymmetricMatrix
    squareform
    """

    def __init__(
        self,
        data: Union[
            np.ndarray,
            Sequence[float],
            Sequence[Sequence[float]],
            "PairwiseMatrix",
        ],
        ids: Optional[Sequence[str]] = None,
        validate: bool = True,
        condensed: bool = False,
        diagonal: Union[float, np.ndarray] = None,
    ):
        data, ids = self._normalize_input(data, ids)

        if ids is None:
            ids = self._generate_ids(data)
        else:
            ids = tuple(ids)

        if validate:
            self._validate_shape(data)
            self._validate_ids(data, ids)
            self._validate_data(data)
            self._validate_diagonal(data, diagonal)

        self._ids = ids
        self._id_index = self._index_list(self._ids)
        self._diagonal = self._init_diagonal(diagonal, data)
        self._data = self._init_data(data, condensed)
        self._flags = self._init_flags(condensed)

    def _init_diagonal(self, diagonal: Union[float, np.ndarray], data: np.ndarray):
        """Initialize the diagonal attribute.

        Parameters
        ----------
        diagonal : float or np.ndarray
            The value or values defined along the diagonal of the matrix.
        data : np.ndarray
            1 or 2-dimensional array of the data of the matrix.

        Returns
        -------
        float or np.ndarray
            The diagonal values of the matrix.
        """
        if diagonal is None:
            if data.ndim == 1:
                return 0.0
            if data.ndim == 2:
                diagonal = np.diagonal(data)
                if np.allclose(diagonal, 0):
                    diagonal = 0.0
                return diagonal
        if np.isscalar(diagonal):
            return float(diagonal)
        else:
            return np.asarray(diagonal)

    def _init_flags(self, condensed: bool) -> dict:
        """Initialize flags for symmetric matrix.

        Parameters
        ----------
        condensed : bool
            Whether the matrix is in condensed form or not.

        Returns
        -------
        dict
            Dictionary containing information on the form of the matrix."""
        if condensed:
            return {"CONDENSED": True}
        else:
            return {"CONDENSED": False}

    def _init_data(self, data: np.ndarray, condensed: bool) -> None:
        """Initialize data for symmetric matrix.

        Parameters
        ----------
        data : np.ndarray
            1 or 2-dimensional array representing the data of the matrix.
        condensed : bool
            Whether or not to store the data in condensed form.

        Returns
        np.ndarray
            1 or 2-dimensional array containing the values of the matrix."""
        if condensed:
            # case where input is 1d and stays 1d
            if data.ndim == 1:
                return data
            # case where input is 2d and is converted to 1d
            else:
                return squareform(data, force="tovector", checks=False)
        else:
            # case where input is 1d and is converted to 2d.
            if data.ndim == 1:
                mat = squareform(data, force="tomatrix", checks=False)
                np.fill_diagonal(mat, self._diagonal)
                return mat
            # case where input is 2d and stays 2d
            else:
                return data

    def _validate_data(self, data: np.ndarray) -> None:
        """Validate the data array.

        Performs a check for symmetry if data is 2D. If data is 1D it is assumed to
        be symmetric.
        """
        if (data.ndim != 1) and (not is_symmetric(data)):
            raise DistanceMatrixError("Data must be symmetric and cannot contain NaNs.")

    def _validate_diagonal(
        self, data: np.ndarray, diagonal: Optional[np.ndarray] = None
    ) -> None:
        """Validate the diagonal of the matrix.

        Checks that the length of the diagonal matches the shape of the matrix, and
        that the diagonal is 1D if passed as an array.

        """
        if diagonal is not None:
            # if it's a single value, it doesn't need to be validated
            if not np.isscalar(diagonal):
                diagonal = np.array(diagonal)
                # if it's a nd.array it needs to be 1d
                if diagonal.ndim != 1:
                    raise SymmetricMatrixError(
                        f"Diagonal must be 1 dimensional if it is an array. Found "
                        f"{diagonal.ndim} dimensions."
                    )
                # it also needs to match the size of the matrix
                length = diagonal.size
                if data.ndim == 2:
                    if length != data.shape[0]:
                        raise SymmetricMatrixError(
                            f"Length of diagonal ({length}) does not match the shape "
                            f"of the matrix {data.shape}."
                        )
                else:
                    shape = _vec_to_shape(data)
                    if length != shape:
                        raise SymmetricMatrixError(
                            f"Length of diagonal ({length}) does not match the shape "
                            f"of the matrix {(shape, shape)}."
                        )

    def _validate_shape(self, data: np.ndarray):
        """Validate the shape of the input data.

        If the input is 2D, it checks that it is at least 1x1, and that it is square.
        If the input is 1D, it checks that it is a correct length to represent a
        symmetric matrix. If input is not 1D or 2D, it raises an error.
        """
        if data.ndim == 2:
            if 0 in data.shape:
                raise SymmetricMatrixError("Data must be at least 1x1 in size.")
            if data.shape[0] != data.shape[1]:
                raise SymmetricMatrixError(
                    "Data must be square (i.e., have the same number of rows and "
                    "columns)."
                )
        # check that the length of the array is a valid length if 1D
        elif data.ndim == 1:
            # pulled from scipy squareform
            s = data.shape
            # Grab the closest value to the square root of the number
            # of elements times 2 to see if the number of elements
            # is indeed a binomial coefficient.
            d = int(np.ceil(np.sqrt(s[0] * 2)))
            # Check that v is of valid dimensions.
            if d * (d - 1) != s[0] * 2:
                raise SymmetricMatrixError(
                    "Incompatible vector size. It must be a binomial "
                    "coefficient n choose 2 for some integer n >= 2."
                )
        else:
            raise SymmetricMatrixError(
                f"Data must be have either 1 or 2 dimensions. "
                f"Found {data.ndim} dimensions."
            )

        if data.dtype not in (np.float32, np.float64):
            raise PairwiseMatrixError("Data must contain only floating point values.")

    @property
    def T(self) -> "SymmetricMatrix":
        r"""Transpose of the matrix.

        If the matrix is in condensed form, a redundant form matrix will be returned.

        See Also
        --------
        transpose

        """
        return self.transpose()

    def transpose(self) -> "SymmetricMatrix":
        r"""Return the transpose of the matrix.

        If the matrix is in condensed form, a redundant form matrix will be returned.

        Notes
        -----
        A deep copy is returned.

        Returns
        -------
        SymmetricMatrix
            Transpose of the matrix. Will be the same type as `self`.

        """
        # Note: Skip validation, since we assume self was already validated
        return self._copy(transpose=True)

    @property
    def diagonal(self) -> Union[float, np.ndarray]:
        """Diagonal value(s) of the matrix.

        If diagonal is a float, this value is repeated along the diagonal of the
        matrix. If diagonal is an array, it represents the full diagonal of the
        matrix."""
        return self._diagonal

    @classonlymethod
    def from_iterable(
        cls,
        iterable: Iterable[Any],
        metric: Callable,
        key: Optional[Any] = None,
        keys: Optional[Iterable[Any]] = None,
        validate: bool = True,
        diagonal=0.0,
        condensed=False,
    ) -> "PairwiseMatrix":
        """Create PairwiseMatrix from an iterable given a metric.

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
        validate : boolean, optional
            If ``True``, all pairwise relationships are computed, including upper
            and lower triangles and the diagonal. If ``False``, `metric` is
            assumed to be hollow and symmetric and only the lower triangle
            (excluding the diagonal) is computed. Pass ``validate=False`` if
            you are sure `metric` is hollow and symmetric for improved
            performance.
        diagonal : float or np.ndarray, optional
            If float, this value will be used to fill the diagonal of the matrix.
            If array, this array will be used to fill the diagonal of the matrix.
            Defaults to 0.0.
        condensed : bool, optional
            If ``False`` (default), the underlying data structure will be the redundant
            form (2D) of the matrix. If ``True``, the underlying data structure will be
            the condensed form (1D) of the matrix.

        Returns
        -------
        PairwiseMatrix
            The `metric` applied to all pairwise elements in the `iterable`.

        Raises
        ------
        ValueError
            If `key` and `keys` are both provided.

        """
        iterable = list(iterable)
        if key is not None and keys is not None:
            raise ValueError("Cannot use both `key` and `keys` at the same time.")

        keys_ = None
        if key is not None:
            keys_ = [resolve_key(e, key) for e in iterable]
        elif keys is not None:
            keys_ = list(keys)

        dm = np.empty((len(iterable),) * 2)
        # this allows for diagonals which do not match the exact shape of the matrix,
        # np.fill_diagonal will just repeat the array to fill. Not sure if this is
        # what we want here.
        np.fill_diagonal(dm, diagonal)
        if validate:
            for i, a in enumerate(iterable):
                for j, b in enumerate(iterable):
                    dm[i, j] = metric(a, b)
        else:
            # This assumes that metric will return a symmetric matrix. That is, that
            # metric(a, b) is the same as metric(b, a)
            for i, a in enumerate(iterable):
                for j, b in enumerate(iterable[:i]):
                    dm[i, j] = dm[j, i] = metric(a, b)
        return cls(dm, keys_, diagonal=diagonal, condensed=condensed)  # type: ignore[operator]

    def __getitem__(
        self, index: Union[str, tuple[str, str], Any]
    ) -> Union[np.ndarray, float]:
        r"""Slice into data by object ID or numpy indexing.

        Extracts data from the matrix by object ID, a pair of
        IDs, or numpy indexing/slicing.

        Parameters
        ----------
        index : str, two-tuple of str, or numpy index
            `index` can be one of the following forms: an ID, a pair of IDs, or
            a numpy index.

            If `index` is a string, it is assumed to be an ID and a
            ``numpy.ndarray`` row vector is returned for the corresponding ID.
            Note that the ID's row of values is returned, *not* its
            column. If the matrix is symmetric, the two will be identical, but
            this makes a difference if the matrix is asymmetric.

            If `index` is a two-tuple of strings, each string is assumed to be
            an ID and the corresponding matrix element is returned that
            represents the value between the two IDs. Note that the
            order of lookup by ID pair matters if the matrix is asymmetric: the
            first ID will be used to look up the row, and the second ID will be
            used to look up the column. Thus, ``dm['a', 'b']`` may not be the
            same as ``dm['b', 'a']`` if the matrix is asymmetric.

            Otherwise, `index` will be passed through to
            ``PairwiseMatrix.data.__getitem__``, allowing for standard
            indexing of a ``numpy.ndarray`` (e.g., slicing).

        Returns
        -------
        ndarray or scalar
            Indexed data, where return type depends on the form of `index` (see
            description of `index` for more details).

        Raises
        ------
        MissingIDError
            If the ID(s) specified in `index` are not in the matrix.

        Notes
        -----
        The lookup based on ID(s) is quick. NumPy indexing (slicing) on condensed
        form matrices will convert them to redundant, roughly doubling their memory
        requirement.

        """
        if isinstance(index, str):
            row_idx = self.index(index)
            if self._flags["CONDENSED"]:
                return _get_row_from_condensed(
                    self._data, row_idx, self.shape[0], self._diagonal
                )
            else:
                return self._data[row_idx]
        elif isinstance(index, tuple) and self._is_id_pair(index):
            i, j = self.index(index[0]), self.index(index[1])
            if self._flags["CONDENSED"]:
                return _get_element_from_condensed(
                    self._data, i, j, self.shape[0], self._diagonal
                )
            else:
                return self._data[i, j]
        else:
            # NumPy index types are numerous and complex, easier to just
            # ignore them in type checking.
            # revert to redundant form to handle numpy style indexing
            if self._flags["CONDENSED"]:
                return self.redundant_form().__getitem__(index)
            else:
                return self._data.__getitem__(index)  # type: ignore[index]

    def redundant_form(self):
        r"""Return an array of values in redundant format.

        Returns
        -------
        ndarray
            Two-dimensional `numpy.ndarray` of values in redundant format.

        See Also
        --------
        condensed_form
        """
        if self._flags["CONDENSED"]:
            mat = squareform(self._data, force="tomatrix", checks=False)
            np.fill_diagonal(mat, self._diagonal)
            return mat
        else:
            return self._data

    def condensed_form(self) -> np.ndarray:
        r"""Return an array of distances in condensed format.

        Returns
        -------
        data : ndarray
            One-dimensional `numpy.ndarray` of values in condensed format.

        Notes
        -----
        Condensed format is described in [1]_.

        The conversion is not a constant-time operation, though it should be
        relatively quick to perform.

        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

        """
        if self._flags["CONDENSED"]:
            return self._data
        else:
            return squareform(self._data, force="tovector", checks=False)

    def permute(
        self,
        condensed: bool = False,
        seed: Optional["SeedLike"] = None,
    ) -> Union["SymmetricMatrix", np.ndarray]:
        r"""Randomly permute both rows and columns in the matrix.

        Randomly permutes the ordering of rows and columns in the matrix. The
        same permutation is applied to both rows and columns in order to
        maintain symmetry and, if applicable, hollowness. Only the rows/columns in the
        matrix are permuted; the IDs are *not* permuted.

        Parameters
        ----------
        condensed : bool, optional
            If ``True``, return the permuted distance matrix in condensed
            format. Otherwise, return the permuted distance matrix as a new
            ``SymmetricMatrix`` instance. Can only be ``True`` if operating on a
            ``DistanceMatrix``.
        seed : int, Generator or RandomState, optional
            A user-provided random seed or random generator instance. See
            :func:`details <skbio.util.get_rng>`.

            .. versionadded:: 0.6.3

        Returns
        -------
        SymmetricMatrix or ndarray
            Permuted values as a new ``DistanceMatrix`` or as a ``ndarray``
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
        rng = get_rng(seed)
        order = rng.permutation(self.shape[0])

        if condensed:
            # if self.__class__ != DistanceMatrix:
            #     raise TypeError("Only distance matrices can return condensed.")
            if self._flags["CONDENSED"]:
                permuted_condensed = distmat_reorder_condensed_py(self._data, order)
            else:
                permuted_condensed = distmat_reorder_condensed(self._data, order)
            return permuted_condensed
        else:
            # Note: Skip validation, since we assume self was already validated
            if self._flags["CONDENSED"]:
                permuted = distmat_reorder_condensed_py(self._data, order)
                return self.__class__(
                    permuted, self.ids, validate=False, condensed=True
                )
            else:
                permuted = distmat_reorder(self._data, order)
                return self.__class__(permuted, self.ids, validate=False)

    def copy(self) -> "SymmetricMatrix":
        r"""Return a deep copy of the symmetric matrix.

        Returns
        -------
        PairwiseMatrix
            Deep copy of the dissimilarity matrix. Will be the same type as
            ``self``.

        """
        if self._flags["CONDENSED"]:
            return self._copy(condensed=True)
        else:
            return self._copy()

    def _copy(
        self, transpose: bool = False, condensed: bool = False
    ) -> "SymmetricMatrix":
        """Copy support.

        Parameters
        ----------
        condensed : bool
            Whether the matrix is in condensed form or not.

        Returns
        -------
        SymmetricMatrix
        """
        # adding for backward compatibility
        data = self._data.copy()
        if transpose:
            data = data.T
        # We deepcopy IDs in case the tuple contains mutable objects at some
        # point in the future.
        # Note: Skip validation, since we assume self was already validated
        return self.__class__(
            data,
            deepcopy(self.ids),
            diagonal=deepcopy(self._diagonal),
            validate=False,
            condensed=condensed,
        )

    def _subset_to_dataframe(
        self, i_ids: Iterable[str], j_ids: Iterable[str]
    ) -> pd.DataFrame:
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
        j_labels = tuple(self.ids[j] for j in j_indices)

        i: list[str] = []
        j: list[str] = []

        # np.hstack([]) throws a ValueError. However, np.hstack([np.array([])])
        # is valid and returns an empty array. Accordingly, an empty array is
        # included here so that np.hstack works in the event that either i_ids
        # or j_ids is empty.
        values = [np.array([])]
        if self._flags["CONDENSED"]:
            diagonal_val = getattr(self, "_diagonal", 0.0)
            n = self.shape[0]

            # precompute to avoid repeated calculations
            condensed_indices = {}
            for i_idx in i_indices:
                for j_idx in j_indices:
                    if i_idx == j_idx:
                        continue
                    key = (min(i_idx, j_idx), max(i_idx, j_idx))
                    if key not in condensed_indices:
                        condensed_indices[key] = _condensed_index(key[0], key[1], n)

            for i_idx in i_indices:
                i.extend([self.ids[i_idx]] * j_length)
                j.extend(j_labels)
                subset_values = np.zeros(j_length, dtype=self._data.dtype)
                for idx, j_idx in enumerate(j_indices):
                    if i_idx == j_idx:
                        if np.isscalar(diagonal_val):
                            subset_values[idx] = diagonal_val
                        else:
                            subset_values[idx] = diagonal_val[i_idx]
                    else:
                        key = (min(i_idx, j_idx), max(i_idx, j_idx))
                        subset_values[idx] = self._data[condensed_indices[key]]
                values.append(subset_values)
        # redundant form
        else:
            for i_idx in i_indices:
                i.extend([self.ids[i_idx]] * j_length)
                j.extend(j_labels)

                subset = self._data[i_idx, j_indices]
                values.append(subset)

        i = pd.Series(i, name="i", dtype=str)
        j = pd.Series(j, name="j", dtype=str)
        values = pd.Series(np.hstack(values), name="value")

        return pd.concat([i, j, values], axis=1)


class DistanceMatrix(SymmetricMatrix):
    """Store distances between objects.

    A `DistanceMatrix` is a `SymmetricMatrix` with the additional
    requirement that the matrix data is hollow. There are additional methods
    made available that take advantage of this hollowness. The
    :func:`~skbio.stats.distance.PairwiseMatrix.plot` method provides
    convenient built-in plotting functionality.

    Parameters
    ----------
    data : array_like or PairwiseMatrix
        Square, hollow, two-dimensional ``numpy.ndarray`` of distances
        (floats), or a structure that can be converted to a ``numpy.ndarray``
        using ``numpy.asarray`` or a one-dimensional vector of distances
        (floats), as defined by `scipy.spatial.distance.squareform <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html>`_.
        Can instead be a `PairwiseMatrix` (or `DistanceMatrix`) instance, in which case
        the instance's data will be used. Data will be converted to a float ``dtype``
        if necessary. A copy will *not* be made if already a ``numpy.ndarray`` with a
        float ``dtype``.
    ids : sequence of str, optional
        Sequence of strings to be used as object IDs. Must match the number of
        rows/cols in `data`. If ``None`` (the default), IDs will be
        monotonically-increasing integers cast as strings, with numbering
        starting from zero, e.g., ``('0', '1', '2', '3', ...)``.
    validate : bool, optional
        If `validate` is ``True`` (the default) and data is not a
        DistanceMatrix object, the input data will be validated.

    See Also
    --------
    PairwiseMatrix
    SymmetricMatrix

    Notes
    -----
    The distances are stored in redundant (square-form) format [1]_. To
    facilitate use with other scientific Python routines (e.g., scipy), the
    distances can be retrieved in condensed (vector-form) format using
    `condensed_form`.

    `DistanceMatrix` only requires that the distances it stores are symmetric and
    hollow. Checks are *not* performed to ensure the other three metric properties
    hold (non-negativity, identity of indiscernibles, and triangle inequality)
    [2]_. Thus, a `DistanceMatrix` instance can store distances that are not
    metric.

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html
    .. [2] http://planetmath.org/metricspace

    """

    # Override here, used in superclass __str__
    _matrix_element_name: ClassVar[str] = "distance"

    def __init__(
        self,
        data: Union[
            np.ndarray,
            Sequence[float],
            Sequence[Sequence[float]],
            "PairwiseMatrix",
        ],
        ids: Optional[Sequence[str]] = None,
        validate: bool = True,
        condensed: bool = False,
        diagonal: Union[float, np.ndarray] = None,
    ):
        data, ids = self._normalize_input(data, ids)

        if ids is None:
            ids = self._generate_ids(data)
        else:
            ids = tuple(ids)

        if validate:
            self._validate_shape(data)
            self._validate_ids(data, ids)
            self._validate_data(data)
            self._validate_diagonal(data, diagonal)

        self._ids = ids
        self._id_index = self._index_list(self._ids)
        self._diagonal = self._init_diagonal()
        self._data = self._init_data(data, condensed)
        self._flags = self._init_flags(condensed)

    def condensed_form(self) -> np.ndarray:
        r"""Return an array of distances in condensed format.

        Returns
        -------
        ndarray
            One-dimensional `numpy.ndarray` of distances in condensed format.

        Notes
        -----
        Condensed format is described in [1]_.

        The conversion is not a constant-time operation, though it should be
        relatively quick to perform.

        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

        """
        if self._flags["CONDENSED"]:
            return self._data
        else:
            return squareform(self._data, force="tovector", checks=False)

    def _validate_data(self, data: np.ndarray) -> None:
        """Validate the data array.

        Performs a check for symmetry and hollowness.
        """
        # if the input data is 1D, we don't need to check for hollowness or symmetry
        if data.ndim == 2:
            data_sym, data_hol = is_symmetric_and_hollow(data)
            if not data_sym:
                raise DistanceMatrixError(
                    "Data must be symmetric and cannot contain NaNs."
                )
            if not data_hol:
                raise DistanceMatrixError(
                    "Data must  be hollow (i.e., the diagonal can only contain zeros)."
                )

    def _validate_diagonal(
        self, data: np.ndarray, diagonal: Union[float, np.ndarray]
    ) -> None:
        """Validate the diagonal of the matrix.

        In addition to the checks performed in the superclass, it also chacks that the
        diagonal only contains zeros.
        """
        if diagonal is not None:
            if np.isscalar(diagonal):
                if diagonal != 0:
                    raise DistanceMatrixError(
                        "The diagonal of a DistanceMatrix may only contain zeros."
                    )
            else:
                diagonal = np.asarray(diagonal)
                if np.any(diagonal != 0):
                    raise DistanceMatrixError(
                        "The diagonal of a DistanceMatrix may only contain zeros."
                    )
        super()._validate_diagonal(data, diagonal)

    def _init_diagonal(self):
        """Initialize the diagonal attribute for a DistanceMatrix.

        A DistanceMatrix must by definition have only 0's on its diagonal, so the most
        efficient storage will always be 0.0."""
        return 0.0

    def to_series(self) -> pd.Series:
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


def randdm(
    num_objects: int,
    ids: Optional[Sequence[str]] = None,
    constructor: Optional[Type[Union["PairwiseMatrix", "DistanceMatrix"]]] = None,
    random_fn: Optional[Union[int, "Generator", Callable]] = None,
) -> "PairwiseMatrix":
    r"""Generate a distance matrix populated with random distances.

    Using the default ``random_fn``, distances are randomly drawn from a uniform
    distribution over ``[0, 1)``.

    Regardless of ``random_fn``, the resulting distance matrix is guaranteed to
    be symmetric and hollow.

    Parameters
    ----------
    num_objects : int
        The number of objects in the resulting distance matrix. For example, if
        ``num_objects`` is 3, a 3x3 distance matrix will be returned.
    ids : sequence of str or None, optional
        A sequence of strings to be used as IDs. ``len(ids)`` must be equal to
        ``num_objects``. If not provided, IDs will be monotonically-increasing
        integers cast as strings (numbering starts at 1). For example,
        ``('1', '2', '3')``.
    constructor : type, optional
        `PairwiseMatrix` or subclass constructor to use when creating the
        random distance matrix. The returned distance matrix will be of this
        type. By default, a `DistanceMatrix` instance will be returned.
    random_fn : int, np.random.Generator or callable, optional
        Function for generating random values. It must accept (n_rows, n_columns) and
        return a 2D array of float-like. Default is the
        :meth:`random <numpy.random.Generator.random>` method of a NumPy random
        generator. If an integer is provided, a random generator will be constructed
        using this number as the seed.

        .. versionchanged:: 0.6.3
            Switched to NumPy's new random generator. Can accept a random seed or
            random generator instance. The function takes one tuple parameter instead
            of two separate parameters.

    Returns
    -------
    PairwiseMatrix
        `PairwiseMatrix` (or subclass) instance of random distances. Type
        depends on ``constructor``.

    See Also
    --------
    numpy.random.Generator.random

    """
    if constructor is None:
        constructor = DistanceMatrix
    if not callable(random_fn):
        random_fn = get_rng(random_fn).random

    assert callable(random_fn)
    data = np.tril(random_fn((num_objects, num_objects)), -1)
    data += data.T

    if not ids:
        ids = tuple(map(str, range(1, num_objects + 1)))

    return constructor(data, ids)


# helper functions for anosim and permanova


def _preprocess_input_sng(
    ids: Sequence,
    sample_size: int,
    grouping: Union[pd.DataFrame, Sequence],
    column: Optional[str],
) -> tuple:
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
            "Grouping vector size must match the number of IDs in the distance matrix."
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


def _preprocess_input(
    distance_matrix: "DistanceMatrix",
    grouping: Union[pd.DataFrame, Sequence],
    column: Optional[str],
) -> tuple:
    """Compute intermediate results not affected by permutations.

    These intermediate results can be computed a single time for efficiency,
    regardless of grouping vector permutations (i.e., when calculating the
    p-value). These intermediate results are used by both ANOSIM and PERMANOVA.

    Also validates and normalizes input (e.g., converting ``DataFrame`` column
    into grouping vector).

    """
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("Input must be a DistanceMatrix.")
    # The if statements are redundant here if I keep the modifications I've made to
    # self.shape, which take advantage of the self._flags dictionary
    # handle redundant form
    # if distance_matrix.data.ndim == 2:
    sample_size = distance_matrix.shape[0]
    # handle condensed form
    # if distance_matrix.data.ndim == 1:
    #     sample_size = _vec_to_shape(distance_matrix.data)

    num_groups, grouping = _preprocess_input_sng(
        distance_matrix.ids, sample_size, grouping, column
    )

    tri_idxs = np.triu_indices(sample_size, k=1)
    distances = distance_matrix.condensed_form()

    return sample_size, num_groups, grouping, tri_idxs, distances


def _df_to_vector(ids: Sequence, df: pd.DataFrame, column: str) -> list:
    """Return a grouping vector from a ``DataFrame`` column.

    Parameters
    ----------
    ids : Sequence
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
            "One or more IDs in the distance matrix are not in the data frame."
        )
    return grouping.tolist()


def _run_monte_carlo_stats(
    test_stat_function: Callable,
    grouping: Union[pd.DataFrame, Sequence],
    permutations: int,
    seed: Optional["SeedLike"] = None,
) -> tuple:
    """Run stat test and compute significance with Monte Carlo permutations."""
    if permutations < 0:
        raise ValueError(
            "Number of permutations must be greater than or equal to zero."
        )

    stat = test_stat_function(grouping)
    rng = get_rng(seed)
    p_value = np.nan
    if permutations > 0:
        perm_stats = np.empty(permutations, dtype=np.float64)

        for i in range(permutations):
            perm_grouping = rng.permutation(grouping)
            perm_stats[i] = test_stat_function(perm_grouping)

        p_value = ((perm_stats >= stat).sum() + 1) / (permutations + 1)

    return stat, p_value


def _build_results(
    method_name: str,
    test_stat_name: str,
    sample_size: int,
    num_groups: int,
    stat: float,
    p_value: float,
    permutations: int,
) -> pd.Series:
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


def _vec_to_shape(vec: np.ndarray) -> int:
    r"""Calculate the redundant shape of a matrix given its condensed vector.

    For an :math:`n \times n` symmetric matrix stored in condensed form
    (upper triangular, excluding diagonal), this function recovers the
    dimension :math:`n` from the number of elements :math:`k` in the
    condensed vector.

    The mathematical relationship is:

    .. math::

        k = \frac{n(n-1)}{2}

    Solving for :math:`n` using the quadratic formula:

    .. math::

        n = \frac{1 + \sqrt{1 + 8k}}{2}

    """
    return int((1 + np.sqrt(1 + 8 * len(vec))) / 2)


def _condensed_index(
    i: Union[int, np.ndarray], j: Union[int, np.ndarray], n: int
) -> Union[int, np.ndarray]:
    """Get indices for condensed form from redundant form indices.

    Parameters
    ----------
    i, j : int or np.ndarray
        Matrix coordinates. Can be scalars or arrays.
    n : int
        Sample size of the square matrix.

    Returns
    -------
    int or np.ndarray
        Index/indices in the condensed form vector.
    """
    is_scalar = np.isscalar(i) and np.isscalar(j)

    i_arr = np.asarray(i)
    j_arr = np.asarray(j)

    swap_mask = i_arr > j_arr
    if np.any(swap_mask):
        i_min = np.where(swap_mask, j_arr, i_arr)
        j_max = np.where(swap_mask, i_arr, j_arr)
    else:
        i_min = i_arr
        j_max = j_arr

    result = i_min * n + j_max - ((i_min + 2) * (i_min + 1)) // 2

    if is_scalar:
        return int(result)

    return result


def distmat_reorder_condensed_py(in_mat, reorder_vec):
    """Pure Python implementation of distmat_reorder for condensed matrices.

    Parameters
    ----------
    in_mat : np.ndarray
        1D condensed form of the matrix.
    reorder_vec : np.ndarray
        1D list of permutation indexes.

    Returns
    -------
    np.ndarray
        Condensed matrix.
    """
    reorder_vec = np.asarray(reorder_vec)
    in_mat = np.asarray(in_mat)
    n_original = _vec_to_shape(in_mat)
    n_filtered = len(reorder_vec)

    # -1 here to exclude the diagonal
    out_size = n_filtered * (n_filtered - 1) // 2
    if out_size == 0:
        return np.array([], dtype=in_mat.dtype)

    i_indices, j_indices = np.triu_indices(n_filtered, k=1)
    old_i = reorder_vec[i_indices]
    old_j = reorder_vec[j_indices]

    old_indices = _condensed_index(old_i, old_j, n_original)

    return in_mat[old_indices]


def _get_element_from_condensed(
    condensed_data: np.ndarray, i: int, j: int, n: int, diagonal: np.ndarray
) -> float:
    """Get a single element from condensed storage.

    Parameters
    ----------
    condensed_data : np.ndarray
        1-dimensional vector form of the matrix.
    i : int
        Row index.
    j : int
        Column index.
    n : int
        Sample size of the square matrix.
    diagonal : float or np.ndarray
        If float, this value is repeated along the diagonal of the matrix.
        If array, this array represents the diagonal of the matrix.

    Returns
    -------
    float
        The value at position (i, j)
    """
    if i == j:
        if np.isscalar(diagonal):
            return diagonal
        else:
            return diagonal[i]

    condensed_idx = _condensed_index(i, j, n)
    return condensed_data[condensed_idx]


def _get_row_from_condensed(
    condensed_data: np.ndarray, row_idx: int, n: int, diagonal: np.ndarray
) -> np.ndarray:
    """Extract a full row from condensed storage.

    Parameters
    ----------
    condensed_data : np.ndarray
        Condensed vector form of the matrix.
    row_idx : int
        Row index of the desired row.
    n : int
        Sample size of the square matrix.
    diagonal : np.ndarray
        Array of values found on the diagonal of the matrix.

    Returns
    -------
    np.ndarray
        The row data for the desired index.
    """
    row = np.empty(n, dtype=condensed_data.dtype)

    # fill in diagonal value
    if np.isscalar(diagonal):
        row[row_idx] = diagonal
    else:
        row[row_idx] = diagonal[row_idx]

    # elements before diagonal: j < row_idx
    j_before = np.arange(row_idx)
    if len(j_before) > 0:
        row[j_before] = condensed_data[_condensed_index(j_before, row_idx, n)]

    # elements after diagonal: j > row_idx
    j_after = np.arange(row_idx + 1, n)
    if len(j_after) > 0:
        row[j_after] = condensed_data[_condensed_index(row_idx, j_after, n)]

    return row


# Alias `DissimilarityMatrix` for backward compatibility
# test whether documentation gets built twice with this
DissimilarityMatrix = PairwiseMatrix
DissimilarityMatrixError = PairwiseMatrixError
