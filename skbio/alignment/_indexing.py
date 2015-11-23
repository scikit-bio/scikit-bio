# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import copy

import numpy as np

class _Indexing(object):
    def __init__(self, instance):
        self._obj = instance

    def __getitem__(self, indexable):
        if type(indexable) is tuple:
            if len(indexable) > 2:
                raise ValueError()
            elif len(indexable) > 1:
                return self._handle_both_axes(*indexable)
            else:
                indexable, = indexable
        return self._slice_on_first_axis(self._obj, indexable)

    def _handle_both_axes(self, seq_index, pos_index):
        if seq_index is Ellipsis or seq_index == slice(None):
            # Only slice second axis
            return self._slice_on_second_axis(self._obj, pos_index)
        else:
            r = self._slice_on_first_axis(self._obj, seq_index)
            if type(r) is not type(self._obj):
                # [1, 1] [1, *]
                return r[pos_index]
            else:
                # [*, 1] [*, *]
                return self._slice_on_second_axis(r, pos_index)

    def _slice_on_second_axis(self, obj, indexable):
        if self.is_scalar(indexable, axis=1):
            # [..., 1]
            return self._get_position(obj, indexable)
        else:
            # [..., *]
            return self._slice_positions(obj, indexable)

    def _slice_on_first_axis(self, obj, indexable):
        if self.is_scalar(indexable, axis=0):
            # [1]
            return self._get_sequence(obj, indexable)
        else:
            # [*]
            return self._slice_sequences(obj, indexable)

    def _get_position(self, obj, indexable):
        return obj._get_position_(indexable)

    def _slice_positions(self, obj, indexable):
        indexable = self._convert_iterable_of_slices(indexable)
        return obj._slice_positions_(indexable)

    def _convert_iterable_of_slices(self, indexable):
        if hasattr(indexable, '__iter__') and not isinstance(indexable,
                                                             np.ndarray):
            indexable = list(indexable)
            if np.asarray(indexable).dtype == object:
                indexable = np.r_[tuple(indexable)]

        return indexable


class TabularMSAILoc(_Indexing):
    def is_scalar(self, indexable, axis):
        return np.isscalar(indexable)

    def _get_sequence(self, obj, indexable):
        return obj._get_sequence_(indexable)

    def _slice_sequences(self, obj, indexable):
        indexable = self._convert_iterable_of_slices(indexable)
        return obj._slice_sequences_(indexable)


class TabularMSALoc(_Indexing):
    def is_scalar(self, indexable, axis):
        """
        Sometimes (MultiIndex!) something that looks like a scalar, isn't
        and vice-versa.

        Consider:

        A 0
          1
          2
        B 0
          1
          2

        'A' looks like a scalar, but isn't.
        ('A', 0) doesn't look like a scalar, but it is.
        """
        index = self._obj.index
        complete_key = False
        partial_key = False
        if axis == 0:
            try:
                if hasattr(index, 'levshape'):
                    if type(indexable) is tuple:
                        complete_key = (len(indexable) == len(index.levshape)
                                        and indexable in index)
                    partial_key = (not complete_key) and indexable in index
            except TypeError:  # Unhashable type, no biggie
                pass
        return (np.isscalar(indexable) and not partial_key) or complete_key

    def _get_sequence(self, obj, indexable):
        return obj._get_sequence_loc_(indexable)

    def _slice_sequences(self, obj, indexable):
        return obj._slice_sequences_loc_(indexable)
