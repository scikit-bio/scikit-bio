# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections.abc

import numpy as np

from skbio.util._decorator import experimental
from ._repr import PairwiseDistancesReprBuilder
from ._util import is_id_pair


class PairwiseDistances(collections.abc.Mapping):
    """Store distances between pairs of IDs.

    ``PairwiseDistances`` is an immutable, unordered, dict-like structure
    (``collections.abc.Mapping``) mapping pairs of IDs to distances. The
    distance between a pair of IDs ("between" distance) must be symmetric, such
    that the distance between IDs A and B is the same as the distance between
    IDs B and A. The distance between an ID and itself ("within" distance) must
    be 0.0 (hollow). It is not necessary to specify "within" distances; they
    are inferred from the input. Likewise, for "between" distances it is not
    necessary to specify both ID pair orientations; only one of A vs. B or
    B vs. A must be supplied, and it does not matter which ordering is used.

    A ``PairwiseDistances`` object can be thought of as a ``DistanceMatrix``
    where not all pairwise distances are defined.
    ``DistanceMatrix.from_pairwise_distances`` can be used to construct a
    ``DistanceMatrix`` from an iterable of ``PairwiseDistances`` objects.
    ``skbio.diversity.beta_diversity`` can produce a ``PairwiseDistances``
    object from a set of ID pairs.

    Parameters
    ----------
    distances : dict
        Dictionary mapping pairs of IDs to distances, where each ID pair is a
        two-element tuple of IDs (strings) mapping to a single distance
        (float).

    See Also
    --------
    skbio.diversity.beta_diversity
    skbio.stats.distance.DistanceMatrix.from_pairwise_distances

    Examples
    --------
    Construct a ``PairwiseDistances`` object from pairs of IDs:

    >>> from skbio.stats.distance import PairwiseDistances
    >>> pwdist = PairwiseDistances({('A', 'B'): 0.5,
    ...                             ('A', 'C'): 0.2,
    ...                             ('D', 'D'): 0.0})
    >>> pwdist
    PairwiseDistances
    --------------------------------------
    Stats:
        ID count: 4
        "Between" distances count: 2
        Total "between" distances count: 6
        Percent complete: 33.33%
    --------------------------------------
    'A': 2 "between" distances
    'B': 1 "between" distance
    'C': 1 "between" distance
    'D': 0 "between" distances

    There are four IDs, `'A'`, `'B'`, `'C'`, and `'D'`. A full
    ``DistanceMatrix`` requires six "between" distances, and this
    ``PairwiseDistances`` object defines two of them (33.33% of the full
    ``DistanceMatrix`` is present). The number of distances between each ID and
    the others are also displayed; `'A'` has two "between" distances
    (A vs. B, A vs. C) while `'B'` and `'C'` each only have a single "between"
    distance (B vs. A and C vs. A, respectively). `'D'` does not have any
    "between" distances defined.

    To get the set of IDs defined on the ``PairwiseDistances`` object:

    >>> pwdist.ids == frozenset({'A', 'B', 'C', 'D'})
    True

    ``PairwiseDistances`` is dict-like, so common dictionary operations may be
    performed, such as checking whether an ID pair exists:

    >>> ('A', 'C') in pwdist
    True
    >>> ('C', 'A') in pwdist
    True
    >>> ('B', 'B') in pwdist
    True
    >>> ('D', 'D') in pwdist
    True
    >>> ('D', 'A') in pwdist
    False
    >>> ('A', 'D') in pwdist
    False

    Note that even though C vs. A was not explicitly provided when constructing
    ``pwdist`` above, both A vs. C and C vs. A are present and map to the same
    distance:

    >>> pwdist['A', 'C']
    0.2
    >>> pwdist['C', 'A']
    0.2

    "Within" distances can also be accessed (they are always 0.0):

    >>> pwdist['B', 'B']
    0.0
    >>> pwdist['D', 'D']
    0.0

    Other dict-like operations are supported, such as length-checking and
    iteration:

    >>> len(pwdist)
    6
    >>> sorted(pwdist)
    [('A', 'A'), ('A', 'B'), ('A', 'C'), ('B', 'B'), ('C', 'C'), ('D', 'D')]

    Note that both the length of and iteration over ``pwdist`` includes both
    "within" and "between" distances.

    """

    # Make object unhashable. May want to revisit this in the future if there
    # is a need.
    __hash__ = None

    # Make object irreversible so that calling reversed(d) raises a TypeError
    # on a Mapping type instead of producing an incorrect iterator. There's an
    # open bug in Python and this is (currently) the recommended way of fixing
    # it: https://bugs.python.org/issue25864
    __reversed__ = None

    @experimental(as_of="0.4.2-dev")
    def __init__(self, distances):
        if not isinstance(distances, dict):
            raise TypeError("`distances` must be a dict, not type %r"
                            % type(distances).__name__)

        distances_ = {}
        ids = set()
        for id_pair, distance in distances.items():
            if not is_id_pair(id_pair):
                raise TypeError(
                    "Each key in `distances` must be a tuple of exactly two "
                    "strings, not %r" % (id_pair,))

            if not isinstance(distance, float):
                raise TypeError(
                    "Each distance in `distances` must be a float, not type %r"
                    % type(distance).__name__)

            if np.isnan(distance):
                raise ValueError("`distances` contains a distance that is NaN")

            id1, id2 = id_pair
            if id1 == id2:
                if distance != 0.0:
                    raise ValueError(
                        "`distances` contains a nonzero distance between ID "
                        "%r and itself: %r" % (id1, distance))
            elif (id2, id1) in distances_:
                if distance != distances_[(id2, id1)]:
                    raise ValueError(
                        "`distances` contains an asymmetric distance between "
                        "IDs %r and %r: %r != %r" % (id1, id2, distance,
                                                     distances_[(id2, id1)]))
            else:
                distances_[id_pair] = distance

            ids.add(id1)
            ids.add(id2)

        self._distances = distances_
        self._ids = frozenset(ids)

    @property
    @experimental(as_of="0.4.2-dev")
    def ids(self):
        """Set of IDs present in this ``PairwiseDistances``.

        Returns
        -------
        frozenset
            Set of IDs present in this ``PairwiseDistances``.

        Notes
        -----
        This property is not writeable.

        Examples
        --------
        >>> from skbio.stats.distance import PairwiseDistances
        >>> pwdist = PairwiseDistances({('A', 'B'): 0.1, ('C', 'C'): 0.0})
        >>> pwdist.ids == frozenset({'A', 'B', 'C'})
        True

        """
        return self._ids

    @experimental(as_of="0.4.2-dev")
    def __bool__(self):
        """Boolean indicating whether this ``PairwiseDistances`` is empty.

        Returns
        -------
        bool
            ``False`` if there are no "within" or "between" distances, ``True``
            otherwise.

        Examples
        --------
        >>> from skbio.stats.distance import PairwiseDistances
        >>> pwdist = PairwiseDistances({('A', 'B'): 0.1, ('C', 'C'): 0.0})
        >>> bool(pwdist)
        True
        >>> pwdist = PairwiseDistances({})
        >>> bool(pwdist)
        False

        """
        return len(self) > 0

    @experimental(as_of="0.4.2-dev")
    def __len__(self):
        """Number of "within" and "between" pairwise distances.

        Returns
        -------
        int
            Number of "within" and "between" pairwise distances in this
            ``PairwiseDistances``.

        Examples
        --------
        >>> from skbio.stats.distance import PairwiseDistances
        >>> pwdist = PairwiseDistances({('A', 'B'): 0.1, ('C', 'C'): 0.0})
        >>> len(pwdist)
        4

        Note that there are three "within" distances and one "between"
        distance.

        """
        return len(self._distances) + len(self.ids)

    @experimental(as_of="0.4.2-dev")
    def __iter__(self):
        """Iterate over "within" and "between" ID pairs in no particular order.

        Yields "within" (ID to same ID) and "between" (ID to different ID)
        ID pairs. Does not yield redundant distances (e.g., only one of A vs.
        B or B vs. A will be yielded; which ordering is yielded is not
        guaranteed).

        There is no guaranteed order to how ID pairs are yielded.

        Yields
        ------
        tuple
            Each "within" and "between" ID pair present in this
            ``PairwiseDistances``.

        Examples
        --------
        >>> from skbio.stats.distance import PairwiseDistances
        >>> pwdist = PairwiseDistances({('A', 'B'): 0.1, ('C', 'C'): 0.0})
        >>> for id_pair in pwdist:
        ...     # do something with `id_pair`, such as `pwdist[id_pair]`
        ...     pass

        """
        yield from self._distances
        for id_ in self.ids:
            yield id_, id_

    @experimental(as_of="0.4.2-dev")
    def __getitem__(self, id_pair):
        """Retrieve the distance between a pair of IDs.

        Supports retrieving both "within" and "between" distances.

        Parameters
        ----------
        id_pair : tuple
            Pair of IDs whose distance will be retrieved. Either ordering of
            IDs is supported, such that A vs. B and B vs. A will return the
            same distance.

        Returns
        -------
        float
            Distance between IDs in `id_pair`.

        Raises
        ------
        TypeError
            If `id_pair` isn't a tuple of exactly two strings.
        KeyError
            If the pair of IDs isn't present in this ``PairwiseDistances``.

        Examples
        --------
        >>> from skbio.stats.distance import PairwiseDistances
        >>> pwdist = PairwiseDistances({('A', 'B'): 0.1, ('C', 'C'): 0.0})
        >>> pwdist['A', 'B']
        0.1
        >>> pwdist['B', 'A']
        0.1
        >>> pwdist['A', 'A']
        0.0
        >>> pwdist['B', 'B']
        0.0
        >>> pwdist['C', 'C']
        0.0

        """
        if is_id_pair(id_pair):
            if id_pair in self._distances:
                return self._distances[id_pair]
            elif (id_pair[1], id_pair[0]) in self._distances:
                return self._distances[id_pair[1], id_pair[0]]
            elif id_pair[0] == id_pair[1] and id_pair[0] in self.ids:
                return 0.0
            else:
                raise KeyError(id_pair)
        else:
            raise TypeError(
                "`id_pair` must be a tuple of exactly two strings, not %r"
                % (id_pair,))

    # A generic implementation of __eq__ is provided by collections.abc.Mapping
    # but we need to provide a different implementation because "between"
    # distances should compare equal in either orientation (A vs. B, B vs. A).
    @experimental(as_of="0.4.2-dev")
    def __eq__(self, other):
        """Determine if this ``PairwiseDistances`` is equal to another.

        ``PairwiseDistances`` objects are equal if they are the same length and
        have the same "within" and "between" pairwise distances.

        Parameters
        ----------
        other : PairwiseDistances
            ``PairwiseDistances`` object to test for equality against.

        Returns
        -------
        bool
            Indicates whether this ``PairwiseDistances`` object is equal to
            `other`.

        See Also
        --------
        __ne__

        Examples
        --------
        >>> from skbio.stats.distance import PairwiseDistances
        >>> pwdist1 = PairwiseDistances({('A', 'B'): 0.1})
        >>> pwdist2 = PairwiseDistances({('B', 'A'): 0.1})
        >>> pwdist1 == pwdist2
        True

        Note that the ordering of ID pairs (A vs. B in ``pwdist1``, B vs. A in
        ``pwdist2``) does not affect the equality of these objects.

        >>> pwdist3 = PairwiseDistances({('B', 'A'): 0.1, ('C', 'C'): 0.0})
        >>> pwdist2 == pwdist3
        False

        """
        if not isinstance(other, PairwiseDistances):
            return False

        if len(self) != len(other):
            return False

        for id_pair in self:
            if id_pair not in other:
                return False
            if self[id_pair] != other[id_pair]:
                return False
        return True

    @experimental(as_of="0.4.2-dev")
    def __ne__(self, other):
        """Determine if this ``PairwiseDistances`` is not equal to another.

        Parameters
        ----------
        other : PairwiseDistances
            ``PairwiseDistances`` object to test for inequality against.

        Returns
        -------
        bool
            Indicates whether this ``PairwiseDistances`` object is not equal to
            `other`.

        See Also
        --------
        __eq__

        """
        return not (self == other)

    @experimental(as_of="0.4.2-dev")
    def __repr__(self):
        """String summary of this ``PairwiseDistances``."""
        pep8_line_length_limit = 79
        length_taken_by_docstring_indent = 8
        width = pep8_line_length_limit - length_taken_by_docstring_indent
        return PairwiseDistancesReprBuilder(
            obj=self,
            width=width,
            indent=4).build()

    def _repr_stats(self):
        id_count = len(self.ids)
        between_count = len(self._distances)
        total_between_count = int((id_count * (id_count - 1)) / 2)
        if id_count > 1:
            percent_complete = between_count / total_between_count
        else:
            # No "between" distances.
            percent_complete = 1.0

        return [
            ("ID count", str(id_count)),
            ('"Between" distances count', str(between_count)),
            ('Total "between" distances count', str(total_between_count)),
            ("Percent complete", '{:.2%}'.format(percent_complete))
        ]

    @experimental(as_of="0.4.2-dev")
    def __str__(self):
        """String summary of this ``PairwiseDistances``."""
        return repr(self)
