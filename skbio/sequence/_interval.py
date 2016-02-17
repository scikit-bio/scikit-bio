# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function
from future.utils import with_metaclass

import abc

from ._feature import Feature
from .intersection import Interval, IntervalTree
import itertools
from collections import defaultdict


class IntervalMetadata():
    def __init__(self, features=None):
        # maps features attributes to intervals
        if features is None:
            self.features = {}
            self.intervals = IntervalTree()
        else:
            self.intervals = IntervalTree()
            for k, invs in features.items():
                for inv in invs:
                    start, end = _polish_interval(inv)
                    self.intervals.add(start, end, k)
            self.features = features

    def add(self, feature, *intervals):
        """ Adds a feature to the metadata object.

        Parameters
        ----------
        feature : skbio.sequence.feature
            The feature object being added.
        intervals : iteerable of intervals
            A list of intervals associated with the feature

        Note
        ----
        The intervals associated with a feature is assumed to be under
        the `location` keyword by default.
        """
        if not isinstance(feature, Feature):
            raise ValueError('feature is not an instance of `Feature`')

        # TODO: The below will require careful consideration
        # since this will be using a little more memory, since the whole
        # feature object will be stored.
        for interval in intervals:
            loc = _polish_interval(interval)
            if loc is not None:
                start, end = loc
                self.intervals.add(start, end, feature)

        # TODO: The below will require careful consideration
        # since this will be using a little more memory, since the whole
        # feature object will be stored.
        if feature not in self.features.keys():
            self.features[feature] = []

        self.features[feature] = list(map(_polish_interval, intervals))


    def _query_interval(self, interval):
        start, end = _polish_interval(interval)
        features = self.intervals.find(start, end)
        return features

    def _query_feature(self, key, value):
        queries = []
        for feature in self.features.keys():
            if feature[key] == value:
                 queries.append(feature)
        return queries

    def query(self, *args, **kwargs):
        """ Looks up features that with intervals and keywords.

        Parameters
        ----------
        args : I1, I2, ...
            Iterable of tuples or Intervals
        kwargs : dict
            Keyword arguments of feature name and feature value, which can
            be passed to ``dict``.  This is used to specify the search parameters.
            If the `location` keyword is passed, then an interval lookup will be
            performed.

        Note
        ----
        There are two types of queries to perform
        1. Query by interval
        2. Query by key/val pair (i.e. gene=sagA)

        Returns
        -------
        list, Feature
            A list of features satisfying the search criteria.
        """
        feats = []

        # Find queries by interval
        for value in args:
            feats += self._query_interval(value)

        # Find queries by feature attribute
        for (key, value) in kwargs.items():
            feats += self._query_feature(key, value)

        return list(set(feats))

    def concat(self, other):
        """ Concatenates two interval metadata objects

        Parameters
        ----------
        other : IntervalMetadata
            An IntervalMetadata object that is being concatenated with
            the current IntervalMetadata object.

        Returns
        -------
        IntervalMetadata
            Concatenated IntervalMetadata object.
        """
        features = {**self.features, **other.features}
        return IntervalMetadata(features)


class IntervalMetadataMixin(with_metaclass(abc.ABCMeta, object)):
    ''' Store metadata corresponding to Features and Intervals.

    Parameters
    ----------
    features : dict of list
        The dictionary of features and intervals where
        the keys are hashble Feature objects and the values
        are a list of intervals associated with each feature.

    There are two underlying data structures namely
    1. An IntervalTree (interval_metadata).
    2. A hash table indexed by a Feature, where each entry contains a
       list of Intervals(feature_metadata).

    Attributes
    ----------
    interval_metadata
    feature_metadata

    '''
    def __init__(self, features=None):
        self._init_(features)

    def _init_(self, features=None):
        self.interval_metadata = IntervalMetadata(features=features)

    def has_interval_metadata(self):
        return self.interval_metadata is not None

def _polish_interval(interval):
    if isinstance(interval, tuple):
        if len(interval) == 0:
            return None
        start, end = interval
        if (len(interval) != 2 or
            ((not isinstance(start, int)) or
             (not isinstance(end, int)))):
            raise ValueError("`start` and `end` aren't correctly specified")
    elif isinstance(interval, Interval):
        start, end = interval.start, interval.end
    elif isinstance(interval, int):
        start, end = interval, interval + 1
    else:
        raise ValueError('The args must be associated with an `Interval` or '
                         'a tuple when querying')
    return start, end
