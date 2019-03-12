# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from abc import ABCMeta, abstractmethod
import collections

import numpy as np
import pandas as pd


class _Indexing(metaclass=ABCMeta):
    def __init__(self, instance, axis=None):
        self._obj = instance
        self._axis = axis

    def __call__(self, axis=None):
        """Set the axis to index on."""
        # verify axis param, discard value
        self._obj._is_sequence_axis(axis)
        return self.__class__(self._obj, axis=axis)

    def __getitem__(self, indexable):
        if self._axis is not None:
            if self._obj._is_sequence_axis(self._axis):
                return self._slice_on_first_axis(self._obj, indexable)
            else:
                return self._slice_on_second_axis(self._obj, indexable)

        if type(indexable) is tuple:
            if len(indexable) > 2:
                raise ValueError("Can only slice on two axes. Tuple is length:"
                                 " %r" % len(indexable))
            elif len(indexable) > 1:
                return self._handle_both_axes(*indexable)
            else:
                indexable, = indexable

        return self._slice_on_first_axis(self._obj, indexable)

    def _handle_both_axes(self, seq_index, pos_index):
        seq_index = self._convert_ellipsis(seq_index)
        pos_index = self._convert_ellipsis(pos_index)

        if not hasattr(seq_index, '__iter__') and seq_index == slice(None):
            # Only slice second axis
            return self._slice_on_second_axis(self._obj, pos_index)
        else:
            r = self._slice_on_first_axis(self._obj, seq_index)
            if type(r) is self._obj.dtype:
                # [1, 1] [1, *]
                return r[pos_index]
            else:
                # [*, 1] [*, *]
                return self._slice_on_second_axis(r, pos_index)

    def _slice_on_second_axis(self, obj, indexable):
        indexable = self._convert_ellipsis(indexable)
        if self.is_scalar(indexable, axis=1):
            # [..., 1]
            return self._get_position(obj, indexable)
        else:
            # [..., *]
            return self._slice_positions(obj, indexable)

    def _slice_on_first_axis(self, obj, indexable):
        indexable = self._convert_ellipsis(indexable)
        if self.is_scalar(indexable, axis=0):
            # [1]
            return self._get_sequence(obj, indexable)
        else:
            # [*]
            return self._slice_sequences(obj, indexable)

    def _convert_ellipsis(self, indexable):
        if indexable is Ellipsis:
            return slice(None)
        return indexable

    @abstractmethod
    def is_scalar(self, indexable, axis):
        raise NotImplementedError

    @abstractmethod
    def _get_sequence(self, obj, indexable):
        raise NotImplementedError

    @abstractmethod
    def _slice_sequences(self, obj, indexable):
        raise NotImplementedError

    def _get_position(self, obj, indexable):
        return obj._get_position_(indexable)

    def _slice_positions(self, obj, indexable):
        indexable = self._assert_bool_vector_right_size(indexable, axis=1)
        indexable = self._convert_iterable_of_slices(indexable)
        return obj._slice_positions_(indexable)

    def _convert_iterable_of_slices(self, indexable):
        # _assert_bool_vector_right_size will have converted the iterable to
        # an ndarray if it wasn't yet.
        if isinstance(indexable, np.ndarray) and indexable.dtype == object:
            indexable = np.r_[tuple(indexable)]

        return indexable

    def _assert_bool_vector_right_size(self, indexable, axis):
        if isinstance(indexable, np.ndarray):
            pass
        elif hasattr(indexable, '__iter__'):
            indexable = np.asarray(list(indexable))
        else:
            return indexable

        if indexable.dtype == bool and len(indexable) != self._obj.shape[axis]:
            raise IndexError("Boolean index's length (%r) does not match the"
                             " axis length (%r)" % (len(indexable),
                                                    self._obj.shape[axis]))

        return indexable


class TabularMSAILoc(_Indexing):
    def is_scalar(self, indexable, axis):
        return np.isscalar(indexable)

    def _get_sequence(self, obj, indexable):
        return obj._get_sequence_iloc_(indexable)

    def _slice_sequences(self, obj, indexable):
        indexable = self._assert_bool_vector_right_size(indexable, axis=0)
        indexable = self._convert_iterable_of_slices(indexable)
        return obj._slice_sequences_iloc_(indexable)


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
        duplicated_key = False
        if not isinstance(indexable, collections.Hashable):
            return False
        if axis == 0 and self._has_fancy_index():
            try:
                if type(indexable) is tuple:
                    complete_key = (len(indexable) == len(index.levshape) and
                                    indexable in index)
                partial_key = not complete_key and indexable in index
            except TypeError:  # Unhashable type, no biggie
                pass
        if index.has_duplicates:
            # for a given Index object index,
            # index[index.duplicated()].unique() is pandas' recommended
            # replacement for index.get_duplicates(), per the pandas 0.23 docs
            duplicated_key = indexable in index[index.duplicated()].unique()
        return (not duplicated_key and
                ((np.isscalar(indexable) and not partial_key) or complete_key))

    def _get_sequence(self, obj, indexable):
        self._assert_tuple_rules(indexable)
        return obj._get_sequence_loc_(indexable)

    def _slice_sequences(self, obj, indexable):
        self._assert_tuple_rules(indexable)
        if (self._has_fancy_index() and
                type(indexable) is not tuple and
                pd.api.types.is_list_like(indexable) and
                len(indexable) > 0):
            if not self.is_scalar(indexable[0], axis=0):
                raise TypeError("A list is used with complete labels, try"
                                " using a tuple to indicate independent"
                                " selections of a `pd.MultiIndex`.")
            # prevents
            # pd.Series.loc[['x', 'b', 'b', 'a']] from being interepereted as
            # pd.Series.loc[[('a', 0), ('b', 1)]] who knows why it does this.
            elif indexable[0] not in self._obj.index:
                raise KeyError(repr(indexable[0]))
            # pandas acts normal if the first element is actually a scalar

        self._assert_bool_vector_right_size(indexable, axis=0)
        return obj._slice_sequences_loc_(indexable)

    def _assert_tuple_rules(self, indexable):
        # pandas is scary in what it will accept sometimes...
        if type(indexable) is tuple:
            if not self._has_fancy_index():
                # prevents unfriendly errors
                raise TypeError("Cannot provide a tuple to the first axis of"
                                " `loc` unless the MSA's `index` is a"
                                " `pd.MultiIndex`.")
            elif self.is_scalar(indexable[0], axis=0):
                # prevents unreasonable results
                # pd.Series.loc[('a', 0), ('b', 1)] would be interpreted as
                # pd.Series.loc[('a', 1)] which is horrifying.
                raise TypeError("A tuple provided to the first axis of `loc`"
                                " represents a selection for each index of a"
                                " `pd.MultiIndex`; it should not contain a"
                                " complete label.")

    def _has_fancy_index(self):
        return hasattr(self._obj.index, 'levshape')
