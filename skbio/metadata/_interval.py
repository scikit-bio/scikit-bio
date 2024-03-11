# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import operator
import copy
import functools

from ._intersection import IntervalTree
from skbio.util._decorator import classonlymethod


class Interval:
    """Stores the bounds and metadata of an interval feature.

    This class stores an interval feature. An interval feature
    is defined as a sub-region of a biological sequence or sequence
    alignment that is a functional entity, e.g., a gene, a riboswitch,
    an exon, etc. It can span a single contiguous region or multiple
    non-contiguous regions (e.g. multiple exons in a transcript, or
    multiple genes in an operon).

    Parameters
    ----------
    interval_metadata : object
        A reference to the ``IntervalMetadata`` object that this
        ``Interval`` object is associated to.
    bounds : iterable of tuple of int
        Tuples representing start and end coordinates. It is *zero-based*
        numbering. It is always inclusive on start bound and exclusive on
        end bound.
    fuzzy : iterable of tuple of bool, optional
        Tuples representing the fuzziness of each bound coordinates.
        If this isn't specified, then the fuzziness of all bound
        coordinates are ``False``. If any of the coordinate fuzziness
        is ``True``, it indicates that the exact bound point of a
        interval feature is unknown. The bound may begin or end at
        some points outside the specified coordinates. This
        accommodates the location format [1]_ of INSDC.
    metadata : dict, optional
        Dictionary of attributes storing information of the feature
        such as "strand", "gene_name", or "product".

    See Also
    --------
    skbio.metadata.IntervalMetadata

    Notes
    -----
    While the construction of an ``Interval`` object automatically add
    itself to its associated ``IntervalMetadata`` object,
    ``IntervalMetadata.add`` is the typical/easier way to
    create and add it to ``IntervalMetadata``.

    References
    ----------
    .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html#3.4.3

    Examples
    --------
    Hypothetically, let's say we have a gene called "genA" with 10 nt
    as shown in the following diagram.  The second row represents the
    two exons (indicated by "=") on this gene:

    ::

        TGGATTCTGC
        -====--==-
        0123456789

    We can create an ``Interval`` object to represent the exons of the gene:

    >>> from skbio.metadata import Interval, IntervalMetadata
    >>> interval_metadata = IntervalMetadata(10)

    Remember the coordinates are inclusive in lower bound and exclusive on
    upper bound:

    >>> gene = Interval(interval_metadata,
    ...                 bounds=[(1, 5), (7, 9)],
    ...                 metadata={'name': 'genA'})
    >>> gene    # doctest: +ELLIPSIS
    Interval(interval_metadata=..., bounds=[(1, 5), (7, 9)], \
fuzzy=[(False, False), (False, False)], metadata={'name': 'genA'})

    """

    def __init__(self, interval_metadata, bounds, fuzzy=None, metadata=None):
        if not isinstance(interval_metadata, IntervalMetadata):
            raise TypeError(
                "You need to provide an IntervalMetadata"
                "object, not %r" % interval_metadata
            )
        # Intervals
        self._interval_metadata = interval_metadata

        self._bounds_fuzzy_setter(bounds, fuzzy)

        # Metadata
        if metadata is None:
            metadata = {}
        self.metadata = metadata

        # add this interval feature to the associated IntervalMetadata
        self._add()

    def _add(self):
        """Add the current ``Interval`` to the IntervalMetadata object."""
        for bound in self.bounds:
            start, end = bound
            self._interval_metadata._interval_tree.add(start, end, self)
        self._interval_metadata._intervals.append(self)

    def __eq__(self, other):
        """Test if this ``Interval`` object is equal to another.

        The equality is performed by checking if the ``metadata``,
        ``bounds`` and ``fuzzy`` are equal. Since the ``bounds``
        and the ``fuzzy`` are sorted, the permutations of them during
        the ``Interval`` construction or assignment won't matter.

        Parameters
        ----------
        other : Interval
            Interval to test for equality against.

        Returns
        -------
        bool
            Indicates if the two objects are equal.

        """
        return (
            (self.metadata == other.metadata)
            and (self.bounds == other.bounds)
            and (self.fuzzy == other.fuzzy)
        )

    def __ne__(self, other):
        """Test if this ``Interval`` object is not equal to another.

        Parameters
        ----------
        other : Interval
            Interval to test for inequality against.

        Returns
        -------
        bool
            Indicates if the two objects are not equal.

        """
        return not (self == other)

    def __repr__(self):
        """Return a string representation of this ``Interval`` object.

        Returns
        -------
        str
            String representation of this ``Interval`` object.

        """
        if self.dropped:
            s = "{}(dropped=True, bounds={!r}, " "fuzzy={!r}, metadata={!r})"
            return s.format(
                self.__class__.__name__, self.bounds, self.fuzzy, self.metadata
            )
        else:
            s = (
                "{}(interval_metadata=<{!r}>, bounds={!r}, "
                "fuzzy={!r}, metadata={!r})"
            )
            return s.format(
                self.__class__.__name__,
                id(self._interval_metadata),
                self.bounds,
                self.fuzzy,
                self.metadata,
            )

    def drop(self):
        """Drop this ``Interval`` object from interval metadata it links to.

        If the ``Interval`` object is dropped, you can still get values of
        ``bounds``, ``fuzzy``, and ``metadata`` attributes, but you
        can not change their values with the setters.

        See Also
        --------
        skbio.metadata.IntervalMetadata.drop

        """
        if not self.dropped:
            self._interval_metadata.drop([self])

    def _bounds_fuzzy_setter(self, bounds=None, fuzzy=None):
        if self.dropped:
            raise RuntimeError(
                "Cannot change `bounds` or `fuzzy` " "on a dropped Interval object."
            )
        # Casts to `list`, validation, sorting, and setting of `bounds`
        # and `fuzzy` happen here.
        if bounds is not None:
            # check iterability
            try:
                # check iterability
                bounds = list(bounds)
            except TypeError:
                raise TypeError(
                    "Cannot give an non-iterable (%r) " "to `bounds`." % bounds
                )

            # check it is not empty
            if not bounds:
                raise ValueError("Cannot give empty `bounds`.")
            # check each contiguous span is in right format
            for bound in bounds:
                _assert_valid_bound(bound)

            spans = len(bounds)
        else:
            spans = len(self.bounds)

        if fuzzy is not None:
            try:
                fuzzy = list(fuzzy)
            except TypeError:
                raise TypeError(
                    "Cannot give a non-iterable (%r) " "to `fuzzy`." % fuzzy
                )

            if len(fuzzy) != spans:
                raise ValueError(
                    "The length of fuzzy must " "be equal to the length of bounds."
                )

            for fuzzy_i in fuzzy:
                _assert_valid_fuzzy(fuzzy_i)

        if bounds is None:
            # `bounds` and `fuzzy` cannot both be omitted.
            if fuzzy is None:
                raise ValueError("Cannot give `None` to both `bounds` " "and `fuzzy`.")
            # If only `fuzzy` is provided, set `self.fuzzy` and don't
            # change `self.bounds`.
            else:
                self._fuzzy = fuzzy
        else:
            # If only `bounds` is provided, reset `self.fuzzy` to
            # all `False`.
            if fuzzy is None:
                bounds.sort()
                self._check_bounds(bounds)
                self._bounds = bounds
                # reset all the fuzzy to False!!
                del self.fuzzy

            # If both `bounds` and `fuzzy` are provided, set
            # `self.bounds` and `self.fuzzy`.
            else:
                bounds, fuzzy = [list(e) for e in zip(*sorted(zip(bounds, fuzzy)))]
                self._check_bounds(bounds)
                self._bounds = bounds
                self._fuzzy = fuzzy

            self._interval_metadata._is_stale_tree = True

    def _check_bounds(self, bounds):
        """Input `bounds` must be sorted."""
        upper_bound = self._interval_metadata.upper_bound
        lower_bound = self._interval_metadata.lower_bound
        if upper_bound is not None and bounds[-1][-1] > upper_bound:
            raise ValueError(
                "Cannot set `bounds` (%r) with coordinate "
                "larger than upper bound (%r)" % (bounds, upper_bound)
            )

        if bounds[0][0] < lower_bound:
            raise ValueError(
                "Cannot set `bounds` (%r) with coordinate "
                "smaller than lower bound (%r)." % (bounds, lower_bound)
            )

    @property
    def fuzzy(self):
        """The openness of each coordinate.

        This indicates that the exact bound of a interval feature
        is unknown. The bound may begin or end at some points outside
        the specified coordinates. This accommodates the bound format [1]_
        of INSDC.

        References
        ----------
        .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html#3.4.3

        """
        return self._fuzzy

    @fuzzy.setter
    def fuzzy(self, value):
        """Set ``fuzzy``.

        The ``value`` should an iterable matching ``self.bounds``.
        """
        self._bounds_fuzzy_setter(fuzzy=value)

    @fuzzy.deleter
    def fuzzy(self):
        """Delete ``fuzzy``.

        This set all fuzzy to be ``False``.
        """
        if self.dropped:
            raise RuntimeError("Cannot change fuzzy on dropped " "Interval object.")
        self._fuzzy = [(False, False)] * len(self.bounds)

    @property
    def bounds(self):
        """The coordinates of the interval feature.

        It should be a list of tuples of int pair. Each tuple stores
        the start and end coordinates of a span of the interval
        feature. The coordinates are *zero-based*. They are inclusive on
        the start and exclusive on the end.
        """
        return self._bounds

    @bounds.setter
    def bounds(self, value):
        """Set ``bounds``.

        WARNING: setting ``bounds`` will reset ``fuzzy`` value to ``False``.
        This is not totally surprising because it is justifiable your old
        ``fuzzy`` don't fit the new bounds.
        """
        self._bounds_fuzzy_setter(bounds=value)

    @property
    def metadata(self):
        """The metadata of the interval feature.

        It stores the metadata (eg. gene name, function, ID, etc.) of
        the interval feature as a ``dict``.
        """
        return self._metadata

    @metadata.setter
    def metadata(self, value):
        if self.dropped:
            raise RuntimeError("Cannot change metadata on dropped " "Interval object.")
        if not isinstance(value, dict):
            raise TypeError("metadata must be a dict, not %r" % value)
        self._metadata = value

    @metadata.deleter
    def metadata(self):
        """Delete metadata.

        This sets metadata to be empty dict.
        """
        if self.dropped:
            raise RuntimeError("Cannot change metadata to dropped " "Interval object.")
        self._metadata = {}

    @property
    def dropped(self):
        """Boolean value indicating if the ``Interval`` object is dropped.

        If it is dropped, it means it is not associated with IntervalMetadata
        object any more.

        Notes
        -----
        This property is not writable.

        See Also
        --------
        skbio.metadata.Interval.drop
        skbio.metadata.IntervalMetadata.drop

        """
        return self._interval_metadata is None


class IntervalMetadata:
    """Stores the interval features.

    ``IntervalMetadata`` object allows storage, modification, and
    querying of interval features covering a region of a single coordinate
    system. For instance, this can be used to store functional annotations
    about genes across a genome. This object is also applied to the sequence
    alignment.

    This object is typically coupled with another object, such as a
    ``Sequence`` object (or its child class), or a ``TabularMSA`` object.

    Parameters
    ----------
    upper_bound : int or None
        Defines the exclusive upper bound of the interval features. No
        coordinate can be greater than it. ``None``
        means that the coordinate space is unbounded.
    copy_from : IntervalMetadata or None, optional
        Create a new object from the input ``IntervalMetadata`` object by
        shallow copying if it is not ``None``. The upper bound of the new
        object will be updated with the ``upper_bound`` parameter specified.

    Notes
    -----
    This class stores coordinates of all feature bounds into a interval
    tree. It allows the speed up of query-by-bound. The building of
    interval tree is deferred until necessary to save computation. It is
    updated from all coordinates only when you need to fetch info from
    the interval tree.

    When you add a method into this class and if you method need to fetch
    info from ``IntervalMetadata._interval_tree``, you should decorate it with
    ``_rebuild_tree``. This decorator will check if the current interval tree
    is stale and will update it if so. Additionally, if your method add,
    delete, or changes the coordinates of any interval features, you should
    set ``self._is_stale_tree`` to ``True`` at the end of your method to
    indicate the interval tree becomes stale.

    See Also
    --------
    skbio.metadata.Interval

    Examples
    --------
    Let's say we have a sequence of length 10 and want to add annotation
    to it. Create an ``IntervalMetadata`` object:

    >>> from skbio.metadata import Interval, IntervalMetadata
    >>> im = IntervalMetadata(10)

    Let's add annotations of 3 genes:

    >>> im.add(bounds=[(3, 9)],
    ...        metadata={'gene': 'sagB'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., bounds=[(3, 9)], \
fuzzy=[(False, False)], metadata={'gene': 'sagB'})
    >>> im.add(bounds=[(3, 7)],
    ...        metadata={'gene': 'sagC'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., bounds=[(3, 7)], \
fuzzy=[(False, False)], metadata={'gene': 'sagC'})
    >>> im.add(bounds=[(1, 2), (4, 7)],
    ...        metadata={'gene': 'sagA'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., bounds=[(1, 2), (4, 7)], \
fuzzy=[(False, False), (False, False)], metadata={'gene': 'sagA'})

    Show the object representation:

    >>> im    # doctest: +ELLIPSIS
    3 interval features
    -------------------
    Interval(interval_metadata=..., bounds=[(3, 9)], \
fuzzy=[(False, False)], metadata={'gene': 'sagB'})
    Interval(interval_metadata=..., bounds=[(3, 7)], \
fuzzy=[(False, False)], metadata={'gene': 'sagC'})
    Interval(interval_metadata=..., bounds=[(1, 2), (4, 7)], \
fuzzy=[(False, False), (False, False)], metadata={'gene': 'sagA'})

    We can sort the genes by their bounds:

    >>> im.sort()
    >>> im    # doctest: +ELLIPSIS
    3 interval features
    -------------------
    Interval(interval_metadata=..., bounds=[(1, 2), (4, 7)], \
fuzzy=[(False, False), (False, False)], metadata={'gene': 'sagA'})
    Interval(interval_metadata=..., bounds=[(3, 7)], \
fuzzy=[(False, False)], metadata={'gene': 'sagC'})
    Interval(interval_metadata=..., bounds=[(3, 9)], \
fuzzy=[(False, False)], metadata={'gene': 'sagB'})

    Query the genes by bound and/or metadata:

    >>> intvls = im.query([(1, 2)], metadata={'gene': 'foo'})
    >>> list(intvls)
    []
    >>> intvls = im.query([(7, 9)])
    >>> list(intvls)  # doctest: +ELLIPSIS
    [Interval(interval_metadata=..., bounds=[(3, 9)], \
fuzzy=[(False, False)], metadata={'gene': 'sagB'})]
    >>> intvls = im.query(metadata={'gene': 'sagA'})
    >>> intvls = list(intvls)
    >>> intvls  # doctest: +ELLIPSIS
    [Interval(interval_metadata=..., bounds=[(1, 2), (4, 7)], \
fuzzy=[(False, False), (False, False)], metadata={'gene': 'sagA'})]

    Drop the gene(s) we get from query:

    >>> im.drop(intvls)
    >>> im.sort()
    >>> im   # doctest: +ELLIPSIS
    2 interval features
    -------------------
    Interval(interval_metadata=..., bounds=[(3, 7)], \
fuzzy=[(False, False)], metadata={'gene': 'sagC'})
    Interval(interval_metadata=..., bounds=[(3, 9)], \
fuzzy=[(False, False)], metadata={'gene': 'sagB'})

    """

    default_write_format = "gff3"

    def __init__(self, upper_bound, copy_from=None):
        self._upper_bound = upper_bound
        if self.upper_bound is not None:
            if self.upper_bound < self.lower_bound:
                raise ValueError(
                    "Cannot set `upper_bound` (%r) "
                    "smaller than `lower_bound` (%r)"
                    % (self.upper_bound, self.lower_bound)
                )
        # List of Interval objects.
        self._intervals = []

        # IntervalTree object to allow faster querying of interval objects.
        self._interval_tree = IntervalTree()

        # Indicates if the IntervalTree needs to be rebuilt.
        self._is_stale_tree = False

        if copy_from is not None:
            for interval in copy_from._intervals:
                bounds_cp = interval.bounds[:]
                fuzzy_cp = interval.fuzzy[:]
                metadata_cp = copy.copy(interval.metadata)

                self.add(bounds_cp, fuzzy=fuzzy_cp, metadata=metadata_cp)

    @property
    def upper_bound(self):
        """The exclusive upper bound of interval features."""
        return self._upper_bound

    @property
    def lower_bound(self):
        """The inclusive lower bound of interval features."""
        return 0

    @property
    def num_interval_features(self):
        """The total number of interval features."""
        return len(self._intervals)

    def _rebuild_tree(method):
        """Rebuild the IntervalTree."""

        @functools.wraps(method)
        def inner(self, *args, **kwargs):
            if self._is_stale_tree is False:
                return method(self, *args, **kwargs)
            self._interval_tree = IntervalTree()
            for f in self._intervals:
                for start, end in f.bounds:
                    self._interval_tree.add(start, end, f)
            self._is_stale_tree = False
            return method(self, *args, **kwargs)

        return inner

    def _reverse(self):
        """Reverse ``IntervalMetadata`` object.

        This operation reverses all of the interval coordinates.
        For instance, this can be used to compare coordinates
        in the forward strand to coordinates in the reversal strand.
        """
        for f in self._intervals:
            try:
                intvls = [
                    (self.upper_bound - x[1], self.upper_bound - x[0])
                    for x in reversed(f.bounds)
                ]
            except TypeError:
                raise TypeError(
                    "You cannot reverse the coordinates "
                    "when the upper bound is `None`"
                )
            f.bounds = intvls

        # DON'T forget this!!!
        self._is_stale_tree = True

    @classonlymethod
    def concat(cls, interval_metadata):
        """Concatenate an iterable of ``IntervalMetadata`` objects.

        It concatenates the multiple ``IntervalMetadata`` objects into
        one coordinate space. The order of the objects in the input
        iterable matters. The coordinate of the second
        ``InterableMetadata`` will be shifted up with the length of
        the first ``IntervalMetadata`` object.

        This function is useful when you concatenate multiple sequences.

        Parameters
        ----------
        interval_metadata : Iterable (IntervalMetadata)
            The interval metadata to concatenate.

        Returns
        -------
        IntervalMetadata
            Concatenated interval metadata.

        Examples
        --------
        >>> from skbio.metadata import IntervalMetadata

        Create two ``IntervalMetadata`` objects:

        >>> im1 = IntervalMetadata(3)
        >>> _ = im1.add([(0, 2)], [(True, False)], {'gene': 'sagA'})
        >>> im2 = IntervalMetadata(4)
        >>> _ = im2.add([(1, 4)], [(True, True)], {'gene': 'sagB'})

        Concatenate them into a single coordinate space. The second
        ``IntervalMetadata``'s interval features are all shifted
        up. The resulting ``IntervalMetadata``'s upper bound is the
        sum of upper bounds of concatenated objects:

        >>> im = IntervalMetadata.concat([im1, im2])
        >>> im   # doctest: +ELLIPSIS
        2 interval features
        -------------------
        Interval(interval_metadata=<...>, bounds=[(0, 2)], \
fuzzy=[(True, False)], metadata={'gene': 'sagA'})
        Interval(interval_metadata=<...>, bounds=[(4, 7)], \
fuzzy=[(True, True)], metadata={'gene': 'sagB'})
        >>> im.upper_bound
        7

        """
        interval_metadata = list(interval_metadata)

        if len(interval_metadata) == 0:
            return cls(0)

        upper_bound = 0
        for im in interval_metadata:
            try:
                upper_bound += im.upper_bound
            except TypeError:
                raise TypeError(
                    "You cannot concat the interval metadata "
                    "because its upper bound is `None`:\n%r" % im
                )
        new = cls(upper_bound)

        length = 0
        for i, im in enumerate(interval_metadata):
            for intvl in im._intervals:
                bounds = intvl.bounds
                fuzzy = intvl.fuzzy
                if i != 0:
                    bounds = [(start + length, end + length) for start, end in bounds]
                new.add(bounds, fuzzy, intvl.metadata)
            length += im.upper_bound

        return new

    def merge(self, other):
        """Merge the interval features of another ``IntervalMetadata`` object.

        It adds all the interval features of the other object into
        ``self``. Note this will not check if there are any duplicates
        of interval features after merge.

        Notes
        -----
        It will raise error if you merge an unbounded
        ``IntervalMetadata`` object to the current object if it is
        bounded. This avoids partially updating the current object if
        the merge fails in the middle of the process due to the
        possibility that some interval features to be merged may live
        outside the current defined upper bound.

        Parameters
        ----------
        other : ``IntervalMetadata``
            The other ``IntervalMetadata`` to be merged.

        """
        if self.upper_bound is not None:
            if other.upper_bound is None:
                raise ValueError(
                    "Cannot merge an unbound IntervalMetadata object "
                    "to a bounded one"
                )
            elif self.upper_bound != other.upper_bound:
                raise ValueError(
                    "The upper bounds of the two IntervalMetadata objects "
                    "are not equal (%r != %r)" % (self.upper_bound, other.upper_bound)
                )
        if self.lower_bound != other.lower_bound:
            raise ValueError(
                "The lower bounds of the two IntervalMetadata objects "
                "are not equal (%d != %d)" % (self.lower_bound, other.lower_bound)
            )
        for intvl in other._intervals:
            self.add(intvl.bounds, intvl.fuzzy, intvl.metadata)

    def sort(self, ascending=True):
        """Sort interval features by their coordinates.

        It sorts by the start coordinate first. If they are the same between
        two interval features, they will be sorted by comparing their end
        coordinates. For example, an interval feature with [(1, 2), (4, 7)]
        will be sorted in front of another one with [(1, 2), (3, 8)].

        Parameters
        ----------
        ascending : bool, optional
            sort in ascending or descending coordinates.

        """
        self._intervals.sort(
            key=lambda i: [i.bounds[0][0], i.bounds[-1][1]], reverse=not ascending
        )

    def add(self, bounds, fuzzy=None, metadata=None):
        """Create and add an ``Interval`` to this ``IntervalMetadata``.

        This method creates an ``Interval`` object and inserts it into
        the ``IntervalMetadata`` object.

        Parameters
        ----------
        bounds : iterable of tuple of ints
            Tuples representing start and end coordinates. It is *zero-based*
            numbering. It is always inclusive on start bound and exclusive on
            end bound.
        fuzzy : iterable of tuple of bool, optional
            Tuples representing the fuzziness of each bound coordinates.
        metadata : dict, optional
            A dictionary of key-value pairs associated with the
            ``Interval`` object.

        Returns
        -------
        Interval
            The ``Interval`` object added.

        See Also
        --------
        skbio.metadata.Interval

        """
        # Add an interval to the tree. Note that the add functionality is
        # built within the Interval constructor.
        return Interval(
            interval_metadata=self, bounds=bounds, fuzzy=fuzzy, metadata=metadata
        )

    @_rebuild_tree
    def _query_interval(self, bound):
        """Yield ``Interval`` objects that overlap with the bound."""
        _assert_valid_bound(bound)

        start, end = bound
        intvls = self._interval_tree.find(start, end)
        # if a ``Interval`` has many non-contiguous spans and
        # multiple of them overlap with the bound, then
        # this ``Interval`` object will be returned
        # multiple times. So we need to remove duplicates.
        seen = set()
        for intvl in intvls:
            if id(intvl) not in seen:
                seen.add(id(intvl))
                yield intvl

    def _query_attribute(self, metadata, intervals=None):
        """Yield ``Interval`` objects based on query attributes.

        Parameters
        ----------
        metadata : dict or ``None``
            If it is ``None``, return empty iterator; if it is
            ``{}``, return an interator of all the ``Interval``
            objects.
        intervals : iterable
            An iterable of ``Interval`` objects.

        """
        if metadata is None:
            return

        if intervals is None:
            intervals = self._intervals

        for intvl in intervals:
            for key, value in metadata.items():
                if key not in intvl.metadata or intvl.metadata[key] != value:
                    break
            else:
                yield intvl

    @_rebuild_tree
    def query(self, bounds=None, metadata=None):
        """Yield ``Interval`` object with the bounds and attributes.

        The ``Interval`` objects must meet both requirements: 1) overlap
        with any of the spans specified by ``bounds``; 2) satisfy
        ``metadata`` specification. For instance, you can identify
        all the recA genes that overlap with (10, 100) or (900, 1000)
        with this code ``interval_metadata.query([(10, 100),
        (900, 1000)], {'gene': 'recA'})``.

        Parameters
        ----------
        bounds : iterable of tuples of int pair, optional
            Specifies bounds to look for the ``Interval``
            objects. An satisfying interval feature only need to overlap with
            one bound. Default (``None``) means all ``Intervals`` meet
            this requirement.

        metadata : dict, optional
            A dictionary of key word attributes associated with the
            ``Interval`` object. It specifies what metadata keywords and
            values to look for. Default (``None``) means all ``Intervals``
            meet this requirement.

        Yields
        ------
        Interval
            ``Interval`` object satisfying the search criteria.

        """
        if bounds is None and metadata is None:
            metadata = {}

        if bounds is None:
            for intvl in self._query_attribute(metadata):
                yield intvl
        else:
            for loc in bounds:
                intvls = self._query_interval(loc)
                if metadata is None:
                    metadata = {}
                for intvl in self._query_attribute(metadata, intvls):
                    yield intvl

    def drop(self, intervals, negate=False):
        """Drop Interval objects.

        The given ``Interval`` objects will be removed and their
        associated ``IntervalMetadata`` will be set to ``None``.

        Parameters
        ----------
        intervals : iterable of ``Interval``
            ``Interval`` objects to drop from this object.
        negate : bool
            Negate the drop operation, i.e. keeping the specified intervals
            instead of dropping them.

        """
        to_delete = {id(f) for f in intervals}

        new_intvls = []
        # iterate through queries and drop them
        for intvl in self._intervals:
            drop = id(intvl) in to_delete
            if negate is True:
                drop = not drop
            if drop:
                intvl._interval_metadata = None
            else:
                new_intvls.append(intvl)

        self._intervals = new_intvls
        self._is_stale_tree = True

    def __eq__(self, other):
        """Test if this object is equal to another.

        It checks if the coordinate spaces are the same between the
        two objects. If so, then check if all the interval features
        are equal between the two objects after sorting them by
        bounds.

        Parameters
        ----------
        other : IntervalMetadata
            Interval metadata to test for equality against.

        Returns
        -------
        bool
            Indicates if the two objects are equal.

        """
        if (
            self.upper_bound != other.upper_bound
            or self.lower_bound != other.lower_bound
        ):
            return False
        else:
            self_intervals = sorted(self._intervals, key=operator.attrgetter("bounds"))
            other_intervals = sorted(
                other._intervals, key=operator.attrgetter("bounds")
            )
            return self_intervals == other_intervals

    def __ne__(self, other):
        """Test if this object is not equal to another.

        Parameters
        ----------
        other : IntervalMetadata
            Interval metadata to test for inequality against.

        Returns
        -------
        bool
            Indicates if the two objects are not equal.

        See Also
        --------
        skbio.metadata.IntervalMetadata.__eq__

        """
        return not (self == other)

    def __repr__(self):
        """Return a string representation of this object.

        Returns
        -------
        str
            String representation of this ``IntervalMetadata`` object.

        """
        n = self.num_interval_features
        l1 = "{} interval feature".format(n)
        if n != 1:
            l1 += "s"
        l2 = "-" * len(l1)

        if n <= 5:
            items = [repr(i) for i in self._intervals]
        else:
            # intentionally overwrite items[2] to make code cleaner
            items = [repr(self._intervals[i]) for i in [0, 1, 2, n - 2, n - 1]]
            items[2] = "..."

        return "\n".join([l1, l2] + items)

    def __copy__(self):
        """Return a shallow copy.

        Notes
        -----
        The ``IntervalMetadata`` copy will have copies of the
        ``Interval`` objects present in this object.  The ``metadata``
        dictionary of each ``Interval`` object will be a shallow copy.

        See Also
        --------
        __deepcopy__

        """
        return self._copy(False, {})

    def __deepcopy__(self, memo):
        """Return a deep copy.

        Notes
        -----
        The ``IntervalMetadata`` copy will have copies of the
        ``Interval`` objects present in this object.  The ``metadata``
        dictionary of each ``Interval`` object will be a deep copy.

        See Also
        --------
        __copy__

        """
        return self._copy(True, memo)

    def _copy(self, deep, memo):
        cp = IntervalMetadata(self.upper_bound)

        for interval in self._intervals:
            # Only need to shallow-copy `bounds` and `fuzzy`
            # because their elements are immutable.
            bounds_cp = interval.bounds[:]
            fuzzy_cp = interval.fuzzy[:]
            if deep:
                metadata_cp = copy.deepcopy(interval.metadata, memo)
            else:
                metadata_cp = copy.copy(interval.metadata)

            cp.add(bounds_cp, fuzzy=fuzzy_cp, metadata=metadata_cp)

        return cp


def _assert_valid_bound(bound):
    if isinstance(bound, tuple):
        try:
            start, end = bound
        except ValueError:
            raise ValueError(
                "A `bound` must be a tuple of exactly "
                "two coordinates, not {!r}".format(bound)
            )
        if not (isinstance(start, int) and isinstance(end, int)) or start > end:
            raise ValueError(
                "`start` (%r) cannot be a larger int " "than `end` (%r)." % (start, end)
            )
    else:
        raise TypeError("Each `bound` must be a tuple, not {!r}".format(bound))


def _assert_valid_fuzzy(fuzzy):
    if isinstance(fuzzy, tuple):
        try:
            start, end = fuzzy
        except ValueError:
            raise ValueError(
                "A `fuzzy` must be a tuple of exactly " "two, not {!r}".format(fuzzy)
            )
        if not (isinstance(start, bool) and isinstance(end, bool)):
            raise TypeError("A `fuzzy` must be a tuple of two booleans")
    else:
        raise TypeError("Each `fuzzy` must be a tuple, not {!r}".format(fuzzy))
