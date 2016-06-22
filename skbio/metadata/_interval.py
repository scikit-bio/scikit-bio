# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._intersection import IntervalTree
from skbio.util._decorator import experimental
from itertools import chain
import operator


class Interval:
    '''Stores the position and metadata of an interval.

    This class supports the storing data for contiguous intervals
    and non-contiguous intervals.

    Parameters
    ----------
    interval_metadata : object
        A reference to the `IntervalMetadata` object that this
        interval is associated to.
    intervals : iterable of tuple of ints
        List of tuples representing start and end coordinates.
    boundaries : iterable of tuple of bool
        List of tuples, representing the openness of each interval.
        If this isn't specified, then all of the boundaries are true.
    metadata : dict
        Dictionary of attributes storing information of the feature
        such as `strand`, `gene_name` or `product`.

    See Also
    --------
    skbio.metadata.IntervalMetadata
    '''
    def __init__(self, interval_metadata, intervals,
                 boundaries=None, metadata=None):

        ivs = list(intervals)

        # Intervals
        self._interval_metadata = interval_metadata
        # http://stackoverflow.com/a/6422754/1167475
        # Used to sort boundaries later

        indices = [i[0] for i in sorted(enumerate(ivs),
                                        key=operator.itemgetter(1))]
        self.intervals = sorted(ivs)

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
            self._metadata = metadata
        else:
            self._metadata = {}

        if interval_metadata is not None:
            self._add_interval()

    def _add_interval(self):
        """ Add an interval object to the IntervalTree within
            IntervalMetadata."""
        # Add directly to the tree.  So no need for _is_stale_tree
        for loc in self.intervals:
            start, end = loc
            self._interval_metadata._intervals.add(start, end, self)
        self._interval_metadata._metadata.append(self)

    # This is required for creating unique sets of intervals
    def _hash(self):
        return hash(tuple(sorted(self.metadata.items()) +
                          self.intervals +
                          self.boundaries))

    def _cmp(self, other):
        """ Comparison operator required for sorting intervals. """
        return self.intervals < other.intervals

    @experimental(as_of='0.4.2-dev')
    def __eq__(self, other):
        return ((self.metadata == other.metadata) and
                (self.intervals == other.intervals) and
                (self.boundaries == other.boundaries))

    @experimental(as_of='0.4.2-dev')
    def __ne__(self, other):
        return not self.__eq__(other)

    @experimental(as_of='0.4.2-dev')
    def __repr__(self):
        # need pformat to print out the dictionary
        # in a consistent manner
        classname = self.__class__.__name__
        num_intervals = len(self.intervals)
        num_keys = len(self.metadata)
        sts, ends = zip(*self.intervals)
        st = min(sts)
        end = max(ends)
        coords = "start=%d, end=%d" % (st, end)

        if num_intervals == 1:
            interval_summary = '1 contiguous interval'
        else:
            interval_summary = '%d non-contiguous intervals' % num_intervals
        summary = ('<%s: %s, %s, %d metadata keys, dropped=%s>' %
                   (classname, coords, interval_summary,
                    num_keys, self.dropped))
        return summary

    @experimental(as_of='0.4.2-dev')
    def __str__(self):
        return self.__repr__()

    @experimental(as_of='0.4.2-dev')
    def drop(self):
        self._interval_metadata.drop(intervals=self.intervals,
                                     boundaries=self.boundaries,
                                     metadata=self.metadata)
        self.boundaries = None
        self.intervals = None
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
    def intervals(self):
        if self.dropped:
            raise RuntimeError('Cannot retrieve intervals from.'
                               'dropped Interval object.')
        return self._intervals

    @intervals.setter
    @experimental(as_of='0.4.2-dev')
    def intervals(self, value):
        if value is None:
            self._intervals = None
        else:
            for interval in value:
                _assert_valid_interval(interval)

            if self.dropped:
                raise RuntimeError('Cannot change intervals to dropped '
                                   'Interval object.')
            self._intervals = value
            if self._interval_metadata is not None:
                self._interval_metadata._is_stale_tree = True

    @property
    @experimental(as_of='0.4.2-dev')
    def metadata(self):
        if self.dropped:
            raise RuntimeError('Cannot retrieve intervals from.'
                               'dropped Interval object.')
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
    """ Stores Interval objects """
    def __init__(self):
        # stores metadata for each feature
        self._metadata = []
        self._intervals = IntervalTree()
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
        for f in self._metadata:
            # staled doesn't need to be called, since the setter for
            # Interval will take care of this
            invs = list(map(lambda x: (length-x[1], length-x[0]),
                            f.intervals))
            f.intervals = invs

    def add(self, intervals, boundaries=None, metadata=None):
        """ Adds a feature to the metadata object.

        Parameters
        ----------
        intervals : iterable of tuple of ints
            A list of intervals associated with the Interval object.
        boundaries : iterable of tuple of bool
            A list of boundaries associated with the Interval object.
        metadata : dict
            A dictionary of key word attributes associated with the
            Interval object.

        Examples
        --------
        >>> from skbio.metadata import IntervalMetadata
        >>> interval_metadata = IntervalMetadata()
        >>> interval_metadata.add(intervals=[(0, 2)],
        ...                       boundaries=None, metadata={'name': 'sagA'})
        >>> interval_metadata.add(intervals=[(40, 70)],
        ...                       boundaries=None, metadata={'name': 'sagB'})
        >>> list(interval_metadata.query(intervals=[(1, 2)]))
        [<Interval: start=0, end=2, 1 contiguous interval, 1 metadata keys, \
dropped=False>]

        """
        # Add an interval to the tree. Note that the add functionality is
        # built within the Interval constructor.
        Interval(interval_metadata=self,
                 intervals=intervals,
                 boundaries=boundaries,
                 metadata=metadata)

    def _rebuild_tree(self, intervals):
        """ Rebuilds the IntervalTree when the tree is stale."""
        self._intervals = IntervalTree()
        for f in intervals:
            for start, end in f.intervals:
                self._intervals.add(start, end, f)

    def _query_interval(self, interval):
        """ Fetches Interval objects based on query interval"""
        _assert_valid_interval(interval)
        start, end = interval
        invs = self._intervals.find(start, end)
        return invs

    def _query_attribute(self, intervals, metadata):
        """ Fetches Interval objects based on query attributes"""
        if metadata is None:
            return []

        for inv in intervals:
            for (key, value) in metadata.items():
                if inv.metadata[key] != value:
                    continue
                yield inv

    @experimental(as_of='0.4.2-dev')
    def query(self, intervals=None, boundaries=None, metadata=None):
        """ Looks up Interval objects with the intervals, boundaries and keywords.

        Parameters
        ----------
        intervals : iterable of tuple of ints
            A list of intervals associated with the Interval object.
        boundaries : iterable of tuple of bool
            A list of boundaries associated with the Interval object.
        metadata : dict
            A dictionary of key word attributes associated with the
            Interval object.

        Returns
        -------
        generator, Interval
            A generator of Interval objects satisfying the search criteria.

        Examples
        --------
        >>> from skbio.metadata import IntervalMetadata
        >>> interval_metadata = IntervalMetadata()
        >>> interval_metadata.add(intervals=[(0, 2)],
        ...                       boundaries=None, metadata={'name': 'sagA'})
        >>> interval_metadata.add(intervals=[(40, 70)],
        ...                       boundaries=None, metadata={'name': 'sagB'})
        >>> list(interval_metadata.query(intervals=[(1, 2)]))
        [<Interval: start=0, end=2, 1 contiguous interval, 1 metadata keys, \
dropped=False>]

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
            self._rebuild_tree(self._metadata)
            self._is_stale_tree = False

        # empty iterator
        def empty():
            return
            yield
        seen = set()
        invs = empty()
        if intervals is None and metadata is None:
            return
            yield
        # only metadata specified
        elif intervals is None and metadata is not None:
            invs = self._query_attribute(self._metadata,
                                         metadata)
            for q in invs:
                if q._hash() not in seen:
                    seen.add(q._hash())
                    yield q

        # only intervals specified
        elif intervals is not None and metadata is None:
            for value in intervals:
                invs = chain(invs, self._query_interval(value))
                for q in invs:
                    if q._hash() not in seen:
                        seen.add(q._hash())
                        yield q
        # both are specified
        else:
            # Find queries by interval
            if intervals is not None:
                for value in intervals:
                    invs = chain(invs, self._query_interval(value))
            invs = self._query_attribute(invs, metadata)
            for q in invs:
                if q._hash() not in seen:
                    seen.add(q._hash())
                    yield q

    @experimental(as_of='0.4.2-dev')
    def drop(self, intervals=None, boundaries=None, metadata=None):
        """ Drops Interval objects according to a specified query.

        Drops all Interval objects that matches the query.

        Parameters
        ----------
        intervals : iterable of tuple of ints
            A list of intervals associated with the Interval object.
        boundaries : iterable of tuple of bool
            A list of boundaries associated with the Interval object.
        metadata : dict
            A dictionary of key word attributes associated with the
            Interval object.

        Examples
        --------
        >>> from skbio.metadata import IntervalMetadata
        >>> interval_metadata = IntervalMetadata()
        >>> interval_metadata.add(intervals=[(0, 2), (4, 7)],
        ...                       boundaries=None, metadata={'name': 'sagA'})
        >>> interval_metadata.add(intervals=[(40, 70)],
        ...                       boundaries=None, metadata={'name': 'sagB'})
        >>> interval_metadata.drop(metadata={'name': 'sagA'})
        >>> list(interval_metadata.query(metadata={'name': 'sagA'}))
        []
        """
        if metadata is None:
            metadata = {}

        queried_invs = self.query(intervals=intervals,
                                  boundaries=boundaries,
                                  metadata=metadata)
        queried_invs = {q._hash(): q for q in queried_invs}

        new_invs = []
        # iterate through queries and drop them
        for inv in self._metadata:
            if inv._hash() not in queried_invs:
                new_invs.append(inv)
        self._metadata = new_invs
        self._is_stale_tree = True

    @experimental(as_of='0.4.2-dev')
    def __eq__(self, other):

        self_metadata = sorted(self._metadata,
                               key=operator.attrgetter('intervals'))
        other_metadata = sorted(other._metadata,
                                key=operator.attrgetter('intervals'))
        return self_metadata == other_metadata

    @experimental(as_of='0.4.2-dev')
    def __ne__(self, other):
        return not self.__eq__(other)

    @experimental(as_of='0.4.2-dev')
    def __repr__(self):
        return str(self._metadata)

    @experimental(as_of='0.4.2-dev')
    def __str__(self):
        return self.__repr__()


def _assert_valid_interval(interval):
    if isinstance(interval, tuple):
        try:
            start, end = interval
            if start > end:
                raise ValueError("`start` is greater than `end`.")
        except:
            raise ValueError("An interval must be a tuple of exactly "
                             "two coordinates, not %r" % (interval, ))
    else:
        raise TypeError("Each interval must be a tuple, not %r" % interval)


def _assert_valid_boundary(boundary):
    if isinstance(boundary, tuple):
        try:
            start, end = boundary
        except:
            raise ValueError("An boundary must be a tuple of exactly "
                             "two coordinates, not %r" % (boundary, ))
    else:
        raise TypeError("Each interval must be a tuple, not %r" % boundary)
