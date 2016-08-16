# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import chain
import operator

from ._intersection import IntervalTree
from skbio.util._decorator import experimental


class Interval:
    """Stores the position and metadata of an interval.

    This class supports the storing data for contiguous intervals
    and non-contiguous intervals.

    Parameters
    ----------
    interval_metadata : object
        A reference to the `IntervalMetadata` object that this
        interval is associated to.
    locations : iterable of tuple of ints
        List of tuples representing start and end coordinates.
    boundaries : iterable of tuple of bool
        List of tuples, representing the openness of each location.
        If this isn't specified, then all of the boundaries are True.
    metadata : dict
        Dictionary of attributes storing information of the feature
        such as `strand`, `gene_name` or `product`.

    Attributes
    ----------
    locations :

    boundaries :

    See Also
    --------
    skbio.metadata.IntervalMetadata


    Examples
    --------
    >>> from skbio.metadata import Interval, IntervalMetadata
    >>> gene1 = Interval(interval_metadata=IntervalMetadata(),
    ...                  locations=[(1, 2), (4, 7)],
    ...                  metadata={'name': 'sagA'})
    >>> gene1    # doctest: +ELLIPSIS
    Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], boundaries=[(True, True), (True, True)], metadata={'name': 'sagA'})
    """
    def __init__(self, interval_metadata, locations,
                 boundaries=None, metadata=None):
        if locations is None:
            raise ValueError('Cannot handle empty set of `locations`.')
        ivs = list(locations)

        # Intervals
        self._interval_metadata = interval_metadata

        # Used to sort boundaries later
        indices = [i[0] for i in sorted(enumerate(ivs),
                                        key=operator.itemgetter(1))]
        self.locations = sorted(ivs)

        # Boundaries
        if boundaries is not None:
            bds = list(boundaries)
            if len(bds) != len(ivs):
                msg = ('The number of boundaries (%d) needs '
                       'to be equal to the number of intervals (%d).')
                raise ValueError(msg % (len(bds), len(ivs)))
            self.boundaries = [bds[i] for i in indices]
        else:
            self.boundaries = [(True, True)] * len(ivs)

        # Metadata
        if metadata is not None:
            self.metadata = metadata
        else:
            self.metadata = {}

        if interval_metadata is not None:
            self._add()

    def _add(self):
        """Add the current `Interval` to the IntervalMetadata object."""
        # Add directly to the tree.  So no need for _is_stale_tree
        # Questions: directly adding leads to rebuilding of the tree?
        for loc in self.locations:
            start, end = loc
            self._interval_metadata._interval_tree.add(start, end, self)
        self._interval_metadata._intervals.append(self)

    def _cmp(self, other):
        """Compare 2 `Interval`s by their locations.

        It is required for sorting by coordinates. Note that this
        method does not check for equality.
        """
        return self.locations < other.locations

    @experimental(as_of='0.4.2-dev')
    def __eq__(self, other):
        ''''''
        return ((self.metadata == other.metadata) and
                (self.locations == other.locations) and
                (self.boundaries == other.boundaries))

    @experimental(as_of='0.4.2-dev')
    def __ne__(self, other):
        ''''''
        return not self.__eq__(other)

    @experimental(as_of='0.4.2-dev')
    def __repr__(self):
        ''''''
        s = '{}(interval_metadata=<{!r}>, locations={!r}, boundaries={!r}, metadata={!r})'
        return s.format(self.__class__.__name__, id(self._interval_metadata),
                        self.locations, self.boundaries, self.metadata)

    @experimental(as_of='0.4.2-dev')
    def __str__(self):
        return self.__repr__()

    @experimental(as_of='0.4.2-dev')
    def drop(self):
        self._interval_metadata.drop(intervals=self.locations,
                                     boundaries=self.boundaries,
                                     metadata=self.metadata)
        self._interval_metadata = None

    @property
    @experimental(as_of='0.4.2-dev')
    def boundaries(self):
        if self.dropped:
            raise RuntimeError('Cannot retrieve boundaries from.'
                               'dropped Interval object.')
        return self._boundaries

    @boundaries.setter
    @experimental(as_of='0.4.2-dev')
    def boundaries(self, value):
        if value is None:
            self._boundaries = None
        else:
            for boundary in value:
                _assert_valid_boundary(boundary)

            if self.dropped:
                raise RuntimeError('Cannot change boundaries to dropped '
                                   'Interval object.')
            self._boundaries = value

    @property
    @experimental(as_of='0.4.2-dev')
    def locations(self):
        if self.dropped:
            raise RuntimeError(
                'Cannot retrieve locations from dropped Interval object.')
        return self._locations

    @locations.setter
    @experimental(as_of='0.4.2-dev')
    def locations(self, value):
        if value is None:
            self._locations = None
        else:
            for location in value:
                _assert_valid_location(location)

            if self.dropped:
                raise RuntimeError(
                    'Cannot change locations to dropped Interval object.')
            self._locations = value
            if self._interval_metadata is not None:
                self._interval_metadata._is_stale_tree = True

    @property
    @experimental(as_of='0.4.2-dev')
    def metadata(self):
        if self.dropped:
            raise RuntimeError(
                'Cannot retrieve metadata from dropped Interval object.')
        return self._metadata

    @metadata.setter
    @experimental(as_of='0.4.2-dev')
    def metadata(self, value):
        if self.dropped:
            raise RuntimeError('Cannot change metadata to dropped '
                               'Interval object.')
        if not isinstance(value, dict):
            raise TypeError("metadata must be a dict, not %r" % value)
        self._metadata = value

    @property
    @experimental(as_of='0.4.2-dev')
    def dropped(self):
        return self._interval_metadata is None


class IntervalMetadata():
    """Stores the interval features of a sequence as a list of `Interval` objects.

    A `IntervalMetadata` includes intervals along a single coordinate
    system. For instance, this can be used to store functional annotations
    about genes across a genome.

    This object is typically coupled with another object, such as a `Sequence`
    object, or a `TabularMSA` object.

    See Also
    --------
    `skbio.metadata.Interval`

    Examples
    --------
    >>> from skbio.metadata import Interval, IntervalMetadata
    >>> from pprint import pprint
    >>> im = IntervalMetadata()
    >>> im.add(locations=[(3, 9)], metadata={'gene': 'sagB'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., locations=[(3, 9)], boundaries=[(True, True)], metadata={'gene': 'sagB'})
    >>> im.add(locations=[(3, 7)], metadata={'gene': 'sagC'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., locations=[(3, 7)], boundaries=[(True, True)], metadata={'gene': 'sagC'})
    >>> im.add(locations=[(1, 2), (4, 7)], metadata={'gene': 'sagA'})  # doctest: +ELLIPSIS
    Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], boundaries=[(True, True), (True, True)], metadata={'gene': 'sagA'})
    >>> im    # doctest: +ELLIPSIS
    3 interval features
    -------------------
    Interval(interval_metadata=..., locations=[(3, 9)], boundaries=[(True, True)], metadata={'gene': 'sagB'})
    Interval(interval_metadata=..., locations=[(3, 7)], boundaries=[(True, True)], metadata={'gene': 'sagC'})
    Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], boundaries=[(True, True), (True, True)], metadata={'gene': 'sagA'})
    >>> im.sort()
    >>> im    # doctest: +ELLIPSIS
    3 interval features
    -------------------
    Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], boundaries=[(True, True), (True, True)], metadata={'gene': 'sagA'})
    Interval(interval_metadata=..., locations=[(3, 7)], boundaries=[(True, True)], metadata={'gene': 'sagC'})
    Interval(interval_metadata=..., locations=[(3, 9)], boundaries=[(True, True)], metadata={'gene': 'sagB'})
    >>> q1 = im.query([(7, 9)])
    >>> list(q1)  # doctest: +ELLIPSIS
    [Interval(interval_metadata=..., locations=[(3, 9)], boundaries=[(True, True)], metadata={'gene': 'sagB'})]
    >>> q2 = im.query([(1, 2), (3, 4)])
    >>> pprint(list(q2))  # doctest: +ELLIPSIS
    [Interval(interval_metadata=..., locations=[(1, 2), (4, 7)], boundaries=[(True, True), (True, True)], metadata={'gene': 'sagA'}),
     Interval(interval_metadata=..., locations=[(3, 7)], boundaries=[(True, True)], metadata={'gene': 'sagC'}),
     Interval(interval_metadata=..., locations=[(3, 9)], boundaries=[(True, True)], metadata={'gene': 'sagB'})]
    >>> q3 = im.query([(1, 2)], metadata={'gene': 'foo'})
    >>> list(q3)
    []
    >>> list(im._query_attributes({'gene': 'foo'}))
    """
    def __init__(self):
        # List of Interval objects.
        self._intervals = []

        # IntervalTree object to allow faster querying of interval objects.
        self._interval_tree = IntervalTree()

        # Indicates if the IntervalTree needs to be rebuilt.
        self._is_stale_tree = False

    def _reverse(self, length):
        """ Reverse complements IntervalMetadata object.

        This operation reverses all of the interval coordinates.
        For instance, this can be used to compare coordinates
        in the forward strand to coordinates in the reversal strand.

        Parameters
        ----------
        length : int
            Largest end coordinate to perform reverse complement.
            This typically corresponds to the length of sequence.
        """
        for f in self._intervals:
            # staled doesn't need to be called, since the setter for
            # Interval will take care of this
            intvls = list(map(lambda x: (length-x[1], length-x[0]),
                              f.locations))
            f.locations = intvls

    def sort(self, ascending=True):
        '''Sort intervals by their starting and ending coordinates.'''
        self._intervals = sorted(
            self._intervals,
            key=lambda i: [i.locations[0][0], i.locations[-1][1]])

    def add(self, locations, boundaries=None, metadata=None):
        """ Adds a feature to the metadata object.

        This method creates a `Interval` object and insert it into
        the `IntervalMetadata` object.

        Parameters
        ----------
        locations : iterable of tuple of ints
            A list of locations associated with the Interval object.
        boundaries : iterable of tuple of bool
            A list of boundaries associated with the locations.
        metadata : dict
            A dictionary of key word attributes associated with the
            Interval object.

        Returns
        -------
        `Interval`
            The `Interval` object just added.

        See Also
        --------
        `skbio.metadata._add`
        """
        # Add an interval to the tree. Note that the add functionality is
        # built within the Interval constructor.
        return Interval(interval_metadata=self,
                        locations=locations,
                        boundaries=boundaries,
                        metadata=metadata)

    def _exist(self, locations, boundaries, metadata):
        '''Return True if the interval data exists.'''
        self._query_interval(locations)

    def _rebuild_tree(self, intervals):
        """ Rebuilds the IntervalTree when the tree is stale."""
        self._interval_tree = IntervalTree()
        for f in intervals:
            for start, end in f.locations:
                self._interval_tree.add(start, end, f)

    def _query_interval(self, location):
        """ Fetches Interval objects that overlap with location."""
        _assert_valid_location(location)
        start, end = location
        intvls = self._interval_tree.find(start, end)
        return intvls

    def _query_attribute(self, metadata, intervals=None):
        """Fetches Interval objects based on query attributes."""
        if metadata is None:
            return []

        if intervals is None:
            intervals = self._intervals

        for intvl in intervals:
            for (key, value) in metadata.items():
                if intvl.metadata[key] != value:
                    break
            else:
                yield intvl

    @experimental(as_of='0.4.2-dev')
    def query(self, locations=None, metadata=None):
        """ Looks up `Interval` objects with the intervals, and keywords.

        All `Interval` objects that satisfy the position constraints and metadata
        will be returned from this function. For instance, this can be used to look for all genes
        within a specific interval in a genome.  Or it could be used to
        find all toxin genes across a subset of a genome.


        Parameters
        ----------
        locations : iterable of tuples of int pair
            A list of intervals associated with the `Interval` object.
            Specifies what ranges of intervals to look for the `Interval`
            objects. An satisfying interval feature only need to overlap

        metadata : dict
            A dictionary of key word attributes associated with the
            Interval object.  Specifies what metadata keywords and
            values to look for.

        Returns
        -------
        generator, Interval
            A generator of Interval objects satisfying the search criteria.

        Note
        ----
        There are two types of queries to perform
        1. Query by interval.
        2. Query by key/val pair (i.e. 'gene':'sagA').

        If no intervals are specified, then only metadata will be searched for.
        If no metadata is specified, then only intervals will be searched for.
        Otherwise, the search will return the intersection of the two results.
        """
        if self._is_stale_tree:
            self._rebuild_tree(self._intervals)
            self._is_stale_tree = False

        # empty iterator
        def empty():
            return
            yield
        seen = set()
        intvls = empty()
        if locations is None and metadata is None:
            return
            yield
        # only metadata specified
        elif locations is None and metadata is not None:
            intvls = self._query_attribute(self._intervals, metadata)
            for q in intvls:
                if id(q) not in seen:
                    seen.add(id(q))
                    yield q

        # only locations specified
        elif locations is not None and metadata is None:
            for value in locations:
                intvls = chain(intvls, self._query_interval(value))
                for q in intvls:
                    if id(q) not in seen:
                        seen.add(id(q))
                        yield q
        # both are specified
        else:
            # Find queries by interval
            if locations is not None:
                for value in locations:
                    intvls = chain(intvls, self._query_interval(value))
            intvls = self._query_attribute(intvls, metadata)
            for q in intvls:
                if id(q) not in seen:
                    seen.add(id(q))
                    yield q

    @experimental(as_of='0.4.2-dev')
    def drop(self, locations=None, boundaries=None, metadata=None):
        """ Drops Interval objects according to a specified query.

        Locations are first queried from the IntervalMetadata object
        using the query functionality. These locations are then dropped
        from the IntervalTree object before being deleted.

        Parameters
        ----------
        locations : iterable of tuple of ints
            A list of locations associated with the Interval object.
        boundaries : iterable of tuple of bool
            A list of boundaries associated with the Interval object.
        metadata : dict
            A dictionary of key word attributes associated with the
            Interval object.

        Examples
        --------
        >>> from skbio.metadata import IntervalMetadata
        >>> im = IntervalMetadata()
        >>> im.add(locations=[(0, 2), (4, 7)], metadata={'name': 'sagA'})   # doctest: +ELLIPSIS
        Interval(interval_metadata=..., locations=[(0, 2), (4, 7)], boundaries=[(True, True), (True, True)], metadata={'name': 'sagA'})
        >>> im.add(locations=[(40, 70)], metadata={'name': 'sagB'})   # doctest: +ELLIPSIS
        Interval(interval_metadata=..., locations=[(40, 70)], boundaries=[(True, True)], metadata={'name': 'sagB'})
        >>> im.drop(metadata={'name': 'sagA'})
        >>> im   # doctest: +ELLIPSIS
        1 interval features
        -------------------
        Interval(interval_metadata=..., locations=[(40, 70)], boundaries=[(True, True)], metadata={'name': 'sagB'})
        """
        if metadata is None:
            metadata = {}

        to_delete = {id(f) for f in
                     self.query(locations=locations,
                                metadata=metadata)}

        new_intvls = []
        # iterate through queries and drop them
        for intvl in self._intervals:
            if id(intvl) not in to_delete:
                new_intvls.append(intvl)

        self._intervals = new_intvls
        self._is_stale_tree = True

    @experimental(as_of='0.4.2-dev')
    def __eq__(self, other):
        ''''''
        self.sort()
        other.sort()
        return self._intervals == other._intervals

    @experimental(as_of='0.4.2-dev')
    def __ne__(self, other):
        return not self.__eq__(other)

    @experimental(as_of='0.4.2-dev')
    def __repr__(self):
        ''''''
        n = len(self._intervals)
        l1 = '{} interval features'.format(n)
        l2 = '-' * len(l1)

        if n <= 5:
            items = [repr(i) for i in self._intervals]
        else:
            items = [repr(self._intervals[i]) for i in [0, 1, n-2, n-1]]
            items[2:2] = '...'

        return '\n'.join([l1, l2, *items])


def _assert_valid_location(location):
    if isinstance(location, tuple):
        try:
            start, end = location
            if start > end:
                raise ValueError("`start` is greater than `end`.")
        except:
            raise ValueError("An location must be a tuple of exactly "
                             "two coordinates, not %r" % (location, ))
    else:
        raise TypeError("Each location must be a tuple, not %r" % location)


def _assert_valid_boundary(boundary):
    if isinstance(boundary, tuple):
        try:
            start, end = boundary
        except:
            raise ValueError("An boundary must be a tuple of exactly "
                             "two coordinates, not %r" % (boundary, ))
    else:
        raise TypeError("Each boundary must be a tuple, not %r" % boundary)
