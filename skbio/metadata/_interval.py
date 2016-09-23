# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import operator
import copy
import functools

from ._intersection import IntervalTree
from skbio.util._decorator import experimental


class Interval:
    """Stores the location and metadata of an interval feature.

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
    locations : iterable of tuple of int
        Tuples representing start and end coordinates.
    boundaries : iterable of tuple of bool, optional
        Tuples representing the openness of each location coordinates.
        If this isn't specified, then all of the boundaries are True. If any
        of the coordinate boundaries is False, it indicates that the exact
        boundary point of a interval feature is unknown. The location may
        begin or end at some points outside the specified coordinates. This
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
    >>> from skbio.metadata import Interval, IntervalMetadata
    >>> interval_metadata = IntervalMetadata(10)

    Create a gene of two exons (from 1 to 2 and 4 to 7):

    >>> gene = Interval(interval_metadata,
    ...                 locations=[(1, 2), (4, 7)],
    ...                 metadata={'name': 'sagA'})
    >>> gene    # doctest: +ELLIPSIS
    Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], \
boundaries=[(True, True), (True, True)], metadata={'name': 'sagA'})
    """
    def __init__(self, interval_metadata, locations,
                 boundaries=None, metadata=None):
        if not isinstance(interval_metadata, IntervalMetadata):
            raise TypeError('You need to provide an IntervalMetadata'
                            'object, not %r' % interval_metadata)
        # Intervals
        self._interval_metadata = interval_metadata

        self._locations_boundaries_setter(locations, boundaries)

        # Metadata
        if metadata is None:
            metadata = {}
        self.metadata = metadata

        # add this interval feature to the associated IntervalMetadata
        self._add()

    def _add(self):
        """Add the current ``Interval`` to the IntervalMetadata object."""
        for loc in self.locations:
            start, end = loc
            self._interval_metadata._interval_tree.add(start, end, self)
        self._interval_metadata._intervals.append(self)

    @experimental(as_of='0.5.0-dev')
    def __eq__(self, other):
        '''Test if this ``Interval`` object is equal to another.

        The equality is performed by checking if the ``metadata``,
        ``location`` and ``boundaries`` are equal. Since the ``locations``
        and the ``boundaries`` are sorted, the permutations of them during
        the ``Interval`` construction or assignment won't matter.

        Parameters
        ----------
        other : Interval
            Interval to test for equality against.

        Returns
        -------
        bool
            Indicates if the two objects are equal.
        '''
        return ((self.metadata == other.metadata) and
                (self.locations == other.locations) and
                (self.boundaries == other.boundaries))

    @experimental(as_of='0.5.0-dev')
    def __ne__(self, other):
        '''Test if this ``Interval`` object is not equal to another.

        Parameters
        ----------
        other : Interval
            Interval to test for inequality against.

        Returns
        -------
        bool
            Indicates if the two objects are not equal.
        '''
        return not self.__eq__(other)

    @experimental(as_of='0.5.0-dev')
    def __repr__(self):
        '''Return a string representation of this ``Interval`` object.

        Returns
        -------
        str
            String representation of this ``Interval`` object.
        '''
        if self.dropped:
            s = ('{}(dropped=True, locations={!r}, '
                 'boundaries={!r}, metadata={!r})')
            return s.format(self.__class__.__name__,
                            self.locations, self.boundaries, self.metadata)
        else:
            s = ('{}(interval_metadata=<{!r}>, locations={!r}, '
                 'boundaries={!r}, metadata={!r})')
            return s.format(self.__class__.__name__,
                            id(self._interval_metadata),
                            self.locations, self.boundaries, self.metadata)

    @experimental(as_of='0.5.0-dev')
    def drop(self):
        '''Drop this ``Interval`` object from the interval metadata it links to.

        If the ``Interval`` object is dropped, you can still get values of
        ``locations``, ``boundaries``, and ``metadata`` attributes, but you
        can not change their values with the setters.

        See Also
        --------
        skbio.metadata.IntervalMetadata.drop
        '''
        if not self.dropped:
            self._interval_metadata.drop([self])

    def _locations_boundaries_setter(self, locations=None, boundaries=None):
        if self.dropped:
            raise RuntimeError('Cannot change `locations` or `boundaries` '
                               'on a dropped Interval object.')
        # Casts to `list`, validation, sorting, and setting of `locations`
        # and `boundaries` happen here.
        if locations is not None:
            # check iterability
            try:
                # check iterability
                locations = list(locations)
            except TypeError:
                raise TypeError('Cannot give an non-iterable (%r) '
                                'to `locations`.' % locations)

            # check it is not empty
            if not locations:
                raise ValueError('Cannot give empty `locations`.')
            # check each contiguous span is in right format
            for location in locations:
                _assert_valid_location(location)

            spans = len(locations)
        else:
            spans = len(self.locations)

        if boundaries is not None:
            try:
                boundaries = list(boundaries)
            except TypeError:
                raise TypeError(
                    'Cannot give a non-iterable (%r) '
                    'to `boundaries`.' % boundaries)

            if len(boundaries) != spans:
                raise ValueError(
                    'The length of boundaries must '
                    'be equal to the length of locations.')

            for boundary in boundaries:
                _assert_valid_boundary(boundary)

        if locations is None:
            # `locations` and `boundaries` cannot both be omitted.
            if boundaries is None:
                raise ValueError('Cannot give `None` to both `locations` '
                                 'and `boundaries`.')
            # If only `boundaries` is provided, set `self.boundaries` and don't
            # change `self.locations`.
            else:
                self._boundaries = boundaries
        else:
            # If only `locations` is provided, reset `self.boundaries` to
            # all `True`.
            if boundaries is None:
                locations.sort()
                self._check_bounds(locations)
                self._locations = locations
                # reset all the boundaries to True!!
                del self.boundaries

            # If both `locations` and `boundaries` are provided, set
            # `self.locations` and `self.boundaries`.
            else:
                locations, boundaries = [
                    list(e) for e in zip(*sorted(zip(locations, boundaries)))]
                self._check_bounds(locations)
                self._locations = locations
                self._boundaries = boundaries

            self._interval_metadata._is_stale_tree = True

    def _check_bounds(self, locations):
        '''`locations must be sorted.'''
        upper_bound = self._interval_metadata.upper_bound
        lower_bound = self._interval_metadata.lower_bound
        if locations[-1][-1] > upper_bound or locations[0][0] < lower_bound:
            raise ValueError('Cannot set `locations` (%r) with coordinate '
                             'larger than upper bound (%r) or smaller than '
                             'lower bound (%r).' %
                             (locations, upper_bound, lower_bound))

    @property
    @experimental(as_of='0.5.0-dev')
    def boundaries(self):
        '''The openness of each coordinate.

        This indicates that the exact boundary point of a interval feature
        is unknown. The location may begin or end at some points outside
        the specified coordinates. This accommodates the location format [1]_
        of INSDC.

        References
        ----------
        .. [1] ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html#3.4.3
        '''
        return self._boundaries

    @boundaries.setter
    @experimental(as_of='0.5.0-dev')
    def boundaries(self, value):
        '''Set ``boundaries``.

        The ``value`` should an iterable matching ``self.locations``.
        '''
        self._locations_boundaries_setter(boundaries=value)

    @boundaries.deleter
    @experimental(as_of='0.5.0-dev')
    def boundaries(self):
        '''Delete ``boundaries``.

        This set all boundaries to be True.
        '''
        if self.dropped:
            raise RuntimeError('Cannot change boundaries on dropped '
                               'Interval object.')
        self._boundaries = [(True, True)] * len(self.locations)

    @property
    @experimental(as_of='0.5.0-dev')
    def locations(self):
        '''The coordinates of the interval feature.

        It should be a list of tuples of int pair. Each tuple stores
        the start and end coordinates of a span of the interval feature.
        '''
        return self._locations

    @locations.setter
    @experimental(as_of='0.5.0-dev')
    def locations(self, value):
        '''Set ``locations``.

        WARNING: setting ``locations`` will reset ``boundaries`` value to True.
        This is not totally surprising because it is justifiable your old
        boundaries don't fit the new locations.
        '''
        self._locations_boundaries_setter(locations=value)

    @property
    @experimental(as_of='0.5.0-dev')
    def metadata(self):
        '''The metadata of the interval feature.

        It stores the metadata (eg. gene name, function, ID, etc.) of
        the interval feature as a ``dict``.
        '''
        return self._metadata

    @metadata.setter
    @experimental(as_of='0.5.0-dev')
    def metadata(self, value):
        if self.dropped:
            raise RuntimeError('Cannot change metadata on dropped '
                               'Interval object.')
        if not isinstance(value, dict):
            raise TypeError("metadata must be a dict, not %r" % value)
        self._metadata = value

    @metadata.deleter
    @experimental(as_of='0.5.0-dev')
    def metadata(self):
        '''Delete metadata.

        This sets metadata to be empty dict.
        '''
        if self.dropped:
            raise RuntimeError('Cannot change metadata to dropped '
                               'Interval object.')
        self._metadata = {}

    @property
    @experimental(as_of='0.5.0-dev')
    def dropped(self):
        '''Boolean value indicating if the ``Interval`` object is dropped.

        If it is dropped, it means it is not associated with IntervalMetadata
        object any more.

        Notes
        -----
        This property is not writable.

        See Also
        --------
        skbio.metadata.Interval.drop
        skbio.metadata.IntervalMetadata.drop
        '''
        return self._interval_metadata is None


class IntervalMetadata():
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
    upper_bound : int
        Defines the upper bound of the interval features. No coordinate can
        be greater than it.

    Notes
    -----
    This class stores coordinates of all feature locations into a interval
    tree. It allows the speed up of query-by-location. The building of
    interval tree is deferred until necessary to save computation. It is
    updated from all coordinates only when you need to fetch info from
    the interval tree.

    When you add a method into this class and if you method need to fetch
    info from ``IntervalMetadata._interval_tree``, you should decorate it with
    ``_rebuild_tree`. This decorator will check if the current interval tree
    is stale and will update it if so. Additionally, if your method add,
    delete, or changes the coordinates of any interval features, you should
    set ``self._is_stale_tree`` to ``True`` at the end of your method to
    indicate the interval tree becomes stale.

    See Also
    --------
    skbio.metadata.Interval

    Examples
    --------
    Create an ``IntervalMetadata`` object:

    >>> from skbio.metadata import Interval, IntervalMetadata
    >>> im = IntervalMetadata(10)

    Let's add some genes annotations:

    >>> im.add(locations=[(3, 9)],
    ...        metadata={'gene': 'sagB'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., locations=[(3, 9)], \
boundaries=[(True, True)], metadata={'gene': 'sagB'})
    >>> im.add(locations=[(3, 7)],
    ...        metadata={'gene': 'sagC'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., locations=[(3, 7)], \
boundaries=[(True, True)], metadata={'gene': 'sagC'})
    >>> im.add(locations=[(1, 2), (4, 7)],
    ...        metadata={'gene': 'sagA'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], \
boundaries=[(True, True), (True, True)], metadata={'gene': 'sagA'})

    Show the object representation:

    >>> im    # doctest: +ELLIPSIS
    3 interval features
    -------------------
    Interval(interval_metadata=..., locations=[(3, 9)], \
boundaries=[(True, True)], metadata={'gene': 'sagB'})
    Interval(interval_metadata=..., locations=[(3, 7)], \
boundaries=[(True, True)], metadata={'gene': 'sagC'})
    Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], \
boundaries=[(True, True), (True, True)], metadata={'gene': 'sagA'})

    We can sort the genes by their locations:

    >>> im.sort()
    >>> im    # doctest: +ELLIPSIS
    3 interval features
    -------------------
    Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], \
boundaries=[(True, True), (True, True)], metadata={'gene': 'sagA'})
    Interval(interval_metadata=..., locations=[(3, 7)], \
boundaries=[(True, True)], metadata={'gene': 'sagC'})
    Interval(interval_metadata=..., locations=[(3, 9)], \
boundaries=[(True, True)], metadata={'gene': 'sagB'})

    Query the genes by location and/or metadata:

    >>> intvls = im.query([(1, 2)], metadata={'gene': 'foo'})
    >>> list(intvls)
    []
    >>> intvls = im.query([(7, 9)])
    >>> list(intvls)  # doctest: +ELLIPSIS
    [Interval(interval_metadata=..., locations=[(3, 9)], \
boundaries=[(True, True)], metadata={'gene': 'sagB'})]
    >>> intvls = im.query(metadata={'gene': 'sagA'})
    >>> intvls = list(intvls)
    >>> intvls  # doctest: +ELLIPSIS
    [Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], \
boundaries=[(True, True), (True, True)], metadata={'gene': 'sagA'})]

    Drop the gene(s) we get from query:

    >>> im.drop(intvls)
    >>> im.sort()
    >>> im   # doctest: +ELLIPSIS
    2 interval features
    -------------------
    Interval(interval_metadata=..., locations=[(3, 7)], \
boundaries=[(True, True)], metadata={'gene': 'sagC'})
    Interval(interval_metadata=..., locations=[(3, 9)], \
boundaries=[(True, True)], metadata={'gene': 'sagB'})
    """
    def __init__(self, upper_bound):
        self._upper_bound = upper_bound
        if self.upper_bound < self.lower_bound:
            raise ValueError('Cannot set `upper_bound` (%r) '
                             'smaller than `lower_bound` (%r)'
                             % (self.upper_bound, self.lower_bound))

        # List of Interval objects.
        self._intervals = []

        # IntervalTree object to allow faster querying of interval objects.
        self._interval_tree = IntervalTree()

        # Indicates if the IntervalTree needs to be rebuilt.
        self._is_stale_tree = False

    @property
    @experimental(as_of='0.5.0-dev')
    def upper_bound(self):
        '''The upper bound of interval features.'''
        return self._upper_bound

    @property
    @experimental(as_of='0.5.0-dev')
    def lower_bound(self):
        '''The lower bound of interval features.'''
        return 0

    @property
    @experimental(as_of='0.5.0-dev')
    def num_interval_features(self):
        '''The total number of interval features.'''
        return len(self._intervals)

    def _rebuild_tree(method):
        """Rebuild the IntervalTree."""
        @functools.wraps(method)
        def inner(self, *args, **kwargs):
            if self._is_stale_tree is False:
                return method(self, *args, **kwargs)
            self._interval_tree = IntervalTree()
            for f in self._intervals:
                for start, end in f.locations:
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
            intvls = [(self.upper_bound - x[1], self.upper_bound - x[0])
                      for x in reversed(f.locations)]
            f.locations = intvls

        # DONT' forget this!!!
        self._is_stale_tree = True

    @experimental(as_of='0.5.0-dev')
    def sort(self, ascending=True):
        '''Sort interval features by their coordinates.

        It sorts by the start coordinate first. If they are the same between
        two interval features, they will be sorted by comparing their end
        coordinates. For example, an interval feature with [(1, 2), (4, 7)]
        will be sorted in front of another one with [(1, 2), (3, 8)].

        Parameters
        ----------
        ascending : bool, optional
            sort in ascending or descending coordinates.
        '''
        self._intervals.sort(
            key=lambda i: [i.locations[0][0], i.locations[-1][1]],
            reverse=not ascending)

    @experimental(as_of='0.5.0-dev')
    def add(self, locations, boundaries=None, metadata=None):
        """Create and add an ``Interval`` to this ``IntervalMetadata``.

        This method creates an ``Interval`` object and inserts it into
        the ``IntervalMetadata`` object.

        Parameters
        ----------
        locations : iterable of tuple of ints
            Tuples representing start and end coordinates.
        boundaries : iterable of tuple of bool, optional
            Tuples representing the openness of each location coordinates.
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
        return Interval(interval_metadata=self,
                        locations=locations,
                        boundaries=boundaries,
                        metadata=metadata)

    @_rebuild_tree
    def _query_interval(self, location):
        """Yield ``Interval`` objects that overlap with the location."""
        _assert_valid_location(location)

        start, end = location
        intvls = self._interval_tree.find(start, end)
        # if a ``Interval`` has many non-contiguous spans and
        # multiple of them overlap with the location, then
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
        intervals : an iterable of ``Interval`` objects
        """
        if metadata is None:
            return

        if intervals is None:
            intervals = self._intervals

        for intvl in intervals:
            for (key, value) in metadata.items():
                if (key not in intvl.metadata or
                        intvl.metadata[key] != value):
                    break
            else:
                yield intvl

    @experimental(as_of='0.5.0-dev')
    @_rebuild_tree
    def query(self, locations=None, metadata=None):
        """Yield ``Interval`` object with the locations and attributes.

        The ``Interval`` objects must meet both requirements: 1) overlap
        with any of the spans specified by ``locations``; 2) satisfy
        ``metadata`` specification. For instance, you can identify
        all the recA genes that overlap with (10, 100) or (900, 1000)
        with this code ``interval_metadata.query([(10, 100),
        (900, 1000)], {'gene': 'recA'})``.

        Parameters
        ----------
        locations : iterable of tuples of int pair, optional
            Specifies locations to look for the ``Interval``
            objects. An satisfying interval feature only need to overlap with
            one location. Default (``None``) means all ``Interval``s meet
            this requirement.

        metadata : dict, optional
            A dictionary of key word attributes associated with the
            ``Interval`` object. It specifies what metadata keywords and
            values to look for. Default (``None``) means all ``Interval``s
            meet this requirement.

        Yields
        ------
        Interval
            ``Interval`` object satisfying the search criteria.
        """
        if locations is None:
            for intvl in self._query_attribute(metadata):
                yield intvl
        else:
            for loc in locations:
                intvls = self._query_interval(loc)
                if metadata is None:
                    metadata = {}
                for intvl in self._query_attribute(metadata, intvls):
                    yield intvl

    @experimental(as_of='0.5.0-dev')
    def drop(self, intervals):
        """Drops Interval objects.

        The given ``Interval`` objects will be removed and their
        associated ``IntervalMetadata`` will be set to ``None``.

        Parameters
        ----------
        intervals : iterable of ``Interval``
            ``Interval`` objects to drop from this object.

        """
        to_delete = {id(f) for f in intervals}

        new_intvls = []
        # iterate through queries and drop them
        for intvl in self._intervals:
            if id(intvl) in to_delete:
                intvl._interval_metadata = None
            else:
                new_intvls.append(intvl)

        self._intervals = new_intvls
        self._is_stale_tree = True

    @experimental(as_of='0.5.0-dev')
    def __eq__(self, other):
        '''Test if this object is equal to another.

        It checks if the coordinate spaces are the same between the
        two objects. If so, then check if all the interval features
        are equal between the two objects after sorting them by
        locations.

        Parameters
        ----------
        other : Interval
            Interval to test for equality against.

        Returns
        -------
        bool
            Indicates if the two objects are equal.

        '''
        if self.upper_bound != other.upper_bound or \
           self.lower_bound != other.lower_bound:
            return False
        else:
            self_intervals = sorted(self._intervals,
                                    key=operator.attrgetter('locations'))
            other_intervals = sorted(other._intervals,
                                     key=operator.attrgetter('locations'))
            return self_intervals == other_intervals

    @experimental(as_of='0.5.0-dev')
    def __ne__(self, other):
        '''Test if this object is not equal to another.

        Parameters
        ----------
        other : Interval
            Interval to test for inequality against.

        Returns
        -------
        bool
            Indicates if the two objects are not equal.

        See Also
        --------
        skbio.metadata.IntervalMetadata.__eq__
        '''
        return not (self == other)

    @experimental(as_of='0.5.0-dev')
    def __repr__(self):
        '''Return a string representation of this object.

        Returns
        -------
        str
            String representation of this ``IntervalMetadata`` object.
        '''
        n = self.num_interval_features
        l1 = '{} interval feature'.format(n)
        if n > 1:
            l1 += 's'
        l2 = '-' * len(l1)

        if n <= 5:
            items = [repr(i) for i in self._intervals]
        else:
            # intentionally overwrite items[2] to make code cleaner
            items = [repr(self._intervals[i]) for i in [0, 1, 2, n-2, n-1]]
            items[2] = '...'

        return '\n'.join([l1, l2] + items)

    @experimental(as_of='0.5.0-dev')
    def __copy__(self):
        '''Return a shallow copy.

        Notes
        -----
        The ``IntervalMetadata`` copy will have copies of the
        ``Interval`` objects present in this object.  The ``metadata``
        dictionary of each ``Interval`` object will be a shallow copy.

        See Also
        --------
        __deepcopy__
        '''
        return self._copy(False, {})

    @experimental(as_of='0.5.0-dev')
    def __deepcopy__(self, memo):
        '''Return a deep copy.

        Notes
        -----
        The ``IntervalMetadata`` copy will have copies of the
        ``Interval`` objects present in this object.  The ``metadata``
        dictionary of each ``Interval`` object will be a deep copy.

        See Also
        --------
        __copy__
        '''
        return self._copy(True, memo)

    def _copy(self, deep, memo):
        cp = IntervalMetadata(self.upper_bound)

        for interval in self._intervals:
            # Only need to shallow-copy `locations` and `boundaries`
            # because their elements are immutable.
            locations_cp = interval.locations[:]
            boundaries_cp = interval.boundaries[:]
            if deep:
                metadata_cp = copy.deepcopy(interval.metadata, memo)
            else:
                metadata_cp = copy.copy(interval.metadata)

            cp.add(locations_cp,
                   boundaries=boundaries_cp,
                   metadata=metadata_cp)

        return cp


def _assert_valid_location(location):
    if isinstance(location, tuple):
        try:
            start, end = location
        except ValueError:
            raise ValueError("A location must be a tuple of exactly "
                             "two coordinates, not {!r}".format(location))
        if not (isinstance(start, int) and
                isinstance(end, int)) or start > end:
            raise ValueError("`start` must be a smaller int than `end`.")
    else:
        raise TypeError("Each location must be a tuple, not {!r}".format(
            location))


def _assert_valid_boundary(boundary):
    if isinstance(boundary, tuple):
        try:
            start, end = boundary
        except ValueError:
            raise ValueError("A boundary must be a tuple of exactly "
                             "two, not {!r}".format(boundary))
        if not (isinstance(start, bool) and isinstance(end, bool)):
            raise TypeError('A boundary must be a tuple of two booleans')
    else:
        raise TypeError("Each boundary must be a tuple, not {!r}".format(
            boundary))
