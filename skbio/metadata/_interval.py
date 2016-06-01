# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function

from ._feature import Feature
from ._intersection import Interval, IntervalTree
from skbio.util._misc import merge_dicts


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

    def reverse_complement(self, length):
        """ Reverse complements IntervalMetadata object.

        Parameters
        ----------
        length : int
            Largest end coordinate to perform reverse complement.
            This typically corresponds to the length of sequence.

        Returns
        -------
        IntervalMetadata
        """
        rvs_features = {}
        for k, v in self.features.items():
            xs = map(_polish_interval, v)
            rvs_features[k] = list(map(lambda x: (length-x[1], length-x[0]),
                                       xs))
        return IntervalMetadata(rvs_features)

    def update(self, old_feature, new_feature):
        """ Replaces a feature in the IntervalMetadata object.

        Finds the original feature and replaces with a new feature.

        Parameters
        ----------
        old_feature : Feature
            Original feature to search and replace in the IntervalMetadata
        new_feature : Feature
            The new feature to replace the original feature.
        """
        ivs = self.features[old_feature]
        self.features[new_feature] = self.features.pop(old_feature)
        for iv in ivs:
            start, end = _polish_interval(iv)
            self.intervals.update(start, end, old_feature, new_feature)

    def add(self, feature, *intervals):
        """ Adds a feature to the metadata object.

        Parameters
        ----------
        feature : skbio.sequence.feature
            The feature object being added.
        intervals : iterable of intervals
            A list of intervals associated with the feature

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
            be passed to ``dict``.  This is used to specify the search
            parameters. If the `location` keyword is passed, then an interval
            lookup will be performed.

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
        feats = set()

        # Find queries by interval
        for value in args:
            feats.update(self._query_interval(value))

        # Find queries by feature attribute
        for (key, value) in kwargs.items():
            feats.update(self._query_feature(key, value))

        return list(feats)

    def concat(self, other, inplace=False):
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

        Notes
        -----
        If the two IntervalMetadata objects contain the same features,
        the features present in other will be used.
        """
        features = merge_dicts(self.features, other.features)
        if inplace:
            self.__init__(features=features)
        else:
            return IntervalMetadata(features)

    def __eq__(self, other):
        # This doesn't look at the interval trees,
        # since the interval trees are strictly built
        # based on the features.
        sivs = list(map(sorted, self.features.values()))
        oivs = list(map(sorted, other.features.values()))

        equalIntervals = sorted(sivs) == sorted(oivs)
        equalFeatures = self.features.keys() == other.features.keys()

        return equalIntervals and equalFeatures


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
