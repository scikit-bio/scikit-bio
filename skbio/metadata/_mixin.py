# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import copy

import pandas as pd

from skbio.metadata import IntervalMetadata


class MetadataMixin(metaclass=abc.ABCMeta):
    @property
    def metadata(self):
        """``dict`` containing metadata which applies to the entire object.

        Notes
        -----
        This property can be set and deleted. When setting new metadata a
        shallow copy of the dictionary is made.

        Examples
        --------
        .. note:: scikit-bio objects with metadata share a common interface for
           accessing and manipulating their metadata. The following examples
           use scikit-bio's ``Sequence`` class to demonstrate metadata
           behavior. These examples apply to all other scikit-bio objects
           storing metadata.

        Create a sequence with metadata:

        >>> from skbio import Sequence
        >>> seq = Sequence('ACGT', metadata={'description': 'seq description',
        ...                                  'id': 'seq-id'})

        Retrieve metadata:

        >>> print(seq.metadata)
        {'description': 'seq description', 'id': 'seq-id'}

        Update metadata:

        >>> seq.metadata['id'] = 'new-id'
        >>> seq.metadata['pubmed'] = 12345
        >>> print(seq.metadata)
        {'description': 'seq description', 'id': 'new-id', 'pubmed': 12345}

        Set metadata:

        >>> seq.metadata = {'abc': 123}
        >>> seq.metadata
        {'abc': 123}

        Delete metadata:

        >>> seq.has_metadata()
        True
        >>> del seq.metadata
        >>> seq.metadata
        {}
        >>> seq.has_metadata()
        False

        """
        if self._metadata is None:
            # Not using setter to avoid copy.
            self._metadata = {}
        return self._metadata

    @metadata.setter
    def metadata(self, metadata):
        if not isinstance(metadata, dict):
            raise TypeError(
                "metadata must be a dict, not type %r" % type(metadata).__name__
            )
        # Shallow copy.
        self._metadata = metadata.copy()

    @metadata.deleter
    def metadata(self):
        self._metadata = None

    @abc.abstractmethod
    def __init__(self, metadata=None):
        raise NotImplementedError

    def _init_(self, metadata=None):
        if metadata is None:
            # Could use deleter but this is less overhead and needs to be fast.
            self._metadata = None
        else:
            # Use setter for validation and copy.
            self.metadata = metadata

    @abc.abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    def _eq_(self, other):
        # We're not simply comparing self.metadata to other.metadata in order
        # to avoid creating "empty" metadata representations on the objects if
        # they don't have metadata.
        if self.has_metadata() and other.has_metadata():
            return self.metadata == other.metadata
        elif not (self.has_metadata() or other.has_metadata()):
            # Both don't have metadata.
            return True
        else:
            # One has metadata while the other does not.
            return False

    @abc.abstractmethod
    def __ne__(self, other):
        raise NotImplementedError

    def _ne_(self, other):
        return not (self == other)

    @abc.abstractmethod
    def __copy__(self):
        raise NotImplementedError

    def _copy_(self):
        if self.has_metadata():
            return self.metadata.copy()
        else:
            return None

    @abc.abstractmethod
    def __deepcopy__(self, memo):
        raise NotImplementedError

    def _deepcopy_(self, memo):
        if self.has_metadata():
            return copy.deepcopy(self.metadata, memo)
        else:
            return None

    def has_metadata(self):
        """Determine if the object has metadata.

        An object has metadata if its ``metadata`` dictionary is not empty
        (i.e., has at least one key-value pair).

        Returns
        -------
        bool
            Indicates whether the object has metadata.

        Examples
        --------
        .. note:: scikit-bio objects with metadata share a common interface for
           accessing and manipulating their metadata. The following examples
           use scikit-bio's ``Sequence`` class to demonstrate metadata
           behavior. These examples apply to all other scikit-bio objects
           storing metadata.

        >>> from skbio import Sequence
        >>> seq = Sequence('ACGT')
        >>> seq.has_metadata()
        False
        >>> seq = Sequence('ACGT', metadata={})
        >>> seq.has_metadata()
        False
        >>> seq = Sequence('ACGT', metadata={'id': 'seq-id'})
        >>> seq.has_metadata()
        True

        """
        return self._metadata is not None and bool(self.metadata)


class PositionalMetadataMixin(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def _positional_metadata_axis_len_(self):
        """Return length of axis that positional metadata applies to.

        Returns
        -------
        int
            Positional metadata axis length.

        """
        raise NotImplementedError

    @property
    def positional_metadata(self):
        """``pd.DataFrame`` containing metadata along an axis.

        Notes
        -----
        This property can be set and deleted. When setting new positional
        metadata, a shallow copy is made and the ``pd.DataFrame`` index is set
        to ``pd.RangeIndex(start=0, stop=axis_len, step=1)``.

        Examples
        --------
        .. note:: scikit-bio objects with positional metadata share a common
           interface for accessing and manipulating their positional metadata.
           The following examples use scikit-bio's ``DNA`` class to demonstrate
           positional metadata behavior. These examples apply to all other
           scikit-bio objects storing positional metadata.

        Create a DNA sequence with positional metadata:

        >>> from skbio import DNA
        >>> seq = DNA(
        ...     'ACGT',
        ...     positional_metadata={'exons': [True, True, False, True],
        ...                          'quality': [3, 3, 20, 11]})
        >>> seq
        DNA
        -----------------------------
        Positional metadata:
            'exons': <dtype: bool>
            'quality': <dtype: int64>
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 50.00%
        -----------------------------
        0 ACGT

        Retrieve positional metadata:

        >>> seq.positional_metadata
           exons  quality
        0   True        3
        1   True        3
        2  False       20
        3   True       11

        Update positional metadata:

        >>> seq.positional_metadata['gaps'] = seq.gaps()
        >>> seq.positional_metadata
           exons  quality   gaps
        0   True        3  False
        1   True        3  False
        2  False       20  False
        3   True       11  False

        Set positional metadata:

        >>> seq.positional_metadata = {'degenerates': seq.degenerates()}
        >>> seq.positional_metadata # doctest: +NORMALIZE_WHITESPACE
          degenerates
        0       False
        1       False
        2       False
        3       False

        Delete positional metadata:

        >>> seq.has_positional_metadata()
        True
        >>> del seq.positional_metadata
        >>> seq.positional_metadata
        Empty DataFrame
        Columns: []
        Index: [0, 1, 2, 3]
        >>> seq.has_positional_metadata()
        False

        """
        if self._positional_metadata is None:
            # Not using setter to avoid copy.
            self._positional_metadata = pd.DataFrame(
                index=self._get_positional_metadata_index()
            )
        return self._positional_metadata

    @positional_metadata.setter
    def positional_metadata(self, positional_metadata):
        try:
            # Pass copy=True to copy underlying data buffer.
            positional_metadata = pd.DataFrame(positional_metadata, copy=True)
        # Different versions of pandas will raise different error types. We
        # don't really care what the type of the error is, just its message, so
        # a blanket Exception will do.
        except Exception as e:
            raise TypeError(
                "Invalid positional metadata. Must be consumable by "
                "`pd.DataFrame` constructor. Original pandas error message: "
                '"%s"' % e
            )

        num_rows = len(positional_metadata.index)
        axis_len = self._positional_metadata_axis_len_()
        if num_rows != axis_len:
            raise ValueError(
                "Number of positional metadata values (%d) must match the "
                "positional metadata axis length (%d)." % (num_rows, axis_len)
            )

        positional_metadata.index = self._get_positional_metadata_index()
        self._positional_metadata = positional_metadata

    @positional_metadata.deleter
    def positional_metadata(self):
        self._positional_metadata = None

    def _get_positional_metadata_index(self):
        """Create a memory-efficient integer index for positional metadata."""
        return pd.RangeIndex(
            start=0, stop=self._positional_metadata_axis_len_(), step=1
        )

    @abc.abstractmethod
    def __init__(self, positional_metadata=None):
        raise NotImplementedError

    def _init_(self, positional_metadata=None):
        if positional_metadata is None:
            # Could use deleter but this is less overhead and needs to be fast.
            self._positional_metadata = None
        else:
            # Use setter for validation and copy.
            self.positional_metadata = positional_metadata

    @abc.abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    def _eq_(self, other):
        # We're not simply comparing self.positional_metadata to
        # other.positional_metadata in order to avoid creating "empty"
        # positional metadata representations on the objects if they don't have
        # positional metadata.
        if self.has_positional_metadata() and other.has_positional_metadata():
            return self.positional_metadata.equals(other.positional_metadata)
        elif not (self.has_positional_metadata() or other.has_positional_metadata()):
            # Both don't have positional metadata.
            return (
                self._positional_metadata_axis_len_()
                == other._positional_metadata_axis_len_()
            )
        else:
            # One has positional metadata while the other does not.
            return False

    @abc.abstractmethod
    def __ne__(self, other):
        raise NotImplementedError

    def _ne_(self, other):
        return not (self == other)

    @abc.abstractmethod
    def __copy__(self):
        raise NotImplementedError

    def _copy_(self):
        if self.has_positional_metadata():
            # deep=True makes a shallow copy of the underlying data buffer.
            return self.positional_metadata.copy(deep=True)
        else:
            return None

    @abc.abstractmethod
    def __deepcopy__(self, memo):
        raise NotImplementedError

    def _deepcopy_(self, memo):
        if self.has_positional_metadata():
            # `copy.deepcopy` no longer recursively copies contents of the
            # DataFrame, so we must handle the deep copy ourselves.
            # Reference: https://github.com/pandas-dev/pandas/issues/17406
            df = self.positional_metadata
            data_cp = copy.deepcopy(df.values.tolist(), memo)
            return pd.DataFrame(
                data_cp,
                index=df.index.copy(deep=True),
                columns=df.columns.copy(deep=True),
                copy=False,
            )
        else:
            return None

    def has_positional_metadata(self):
        """Determine if the object has positional metadata.

        An object has positional metadata if its ``positional_metadata``
        ``pd.DataFrame`` has at least one column.

        Returns
        -------
        bool
            Indicates whether the object has positional metadata.

        Examples
        --------
        .. note:: scikit-bio objects with positional metadata share a common
           interface for accessing and manipulating their positional metadata.
           The following examples use scikit-bio's ``DNA`` class to demonstrate
           positional metadata behavior. These examples apply to all other
           scikit-bio objects storing positional metadata.

        >>> import pandas as pd
        >>> from skbio import DNA
        >>> seq = DNA('ACGT')
        >>> seq.has_positional_metadata()
        False
        >>> seq = DNA('ACGT', positional_metadata=pd.DataFrame(index=range(4)))
        >>> seq.has_positional_metadata()
        False
        >>> seq = DNA('ACGT', positional_metadata={'quality': range(4)})
        >>> seq.has_positional_metadata()
        True

        """
        return (
            self._positional_metadata is not None
            and len(self.positional_metadata.columns) > 0
        )


class IntervalMetadataMixin(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def _interval_metadata_axis_len_(self):
        """Return length of axis that interval metadata applies to.

        Returns
        -------
        int
            Interval metadata axis length.

        """
        raise NotImplementedError

    @abc.abstractmethod
    def __init__(self, interval_metadata=None):
        raise NotImplementedError

    def _init_(self, interval_metadata=None):
        if interval_metadata is None:
            # Could use deleter but this is less overhead and needs to be fast.
            self._interval_metadata = None
        else:
            # Use setter for validation and copy.
            self.interval_metadata = interval_metadata

    @property
    def interval_metadata(self):
        """``IntervalMetadata`` object containing info about interval features.

        Notes
        -----
        This property can be set and deleted. When setting new
        interval metadata, a shallow copy of the ``IntervalMetadata``
        object is made.

        """
        if self._interval_metadata is None:
            # Not using setter to avoid copy.
            self._interval_metadata = IntervalMetadata(
                self._interval_metadata_axis_len_()
            )
        return self._interval_metadata

    @interval_metadata.setter
    def interval_metadata(self, interval_metadata):
        if isinstance(interval_metadata, IntervalMetadata):
            upper_bound = interval_metadata.upper_bound
            lower_bound = interval_metadata.lower_bound
            axis_len = self._interval_metadata_axis_len_()
            if lower_bound != 0:
                raise ValueError(
                    "The lower bound for the interval features (%d) "
                    "must be zero." % lower_bound
                )
            if upper_bound is not None and upper_bound != axis_len:
                raise ValueError(
                    "The upper bound for the interval features (%d) "
                    "must match the interval metadata axis length (%d)"
                    % (upper_bound, axis_len)
                )
            # copy all the data to the mixin
            self._interval_metadata = IntervalMetadata(
                axis_len, copy_from=interval_metadata
            )
        else:
            raise TypeError(
                "You must provide `IntervalMetadata` object, "
                "not type %s." % type(interval_metadata).__name__
            )

    @interval_metadata.deleter
    def interval_metadata(self):
        self._interval_metadata = None

    def has_interval_metadata(self):
        """Determine if the object has interval metadata.

        An object has interval metadata if its ``interval_metadata``
        has at least one ```Interval`` objects.

        Returns
        -------
        bool
            Indicates whether the object has interval metadata.

        """
        return (
            self._interval_metadata is not None
            and self.interval_metadata.num_interval_features > 0
        )

    @abc.abstractmethod
    def __eq__(self, other):
        raise NotImplementedError

    def _eq_(self, other):
        # We're not simply comparing self.interval_metadata to
        # other.interval_metadata in order to avoid creating "empty"
        # interval metadata representations on the objects if they don't have
        # interval metadata.
        if self.has_interval_metadata() and other.has_interval_metadata():
            return self.interval_metadata == other.interval_metadata
        elif not (self.has_interval_metadata() or other.has_interval_metadata()):
            # Both don't have interval metadata.
            return (
                self._interval_metadata_axis_len_()
                == other._interval_metadata_axis_len_()
            )
        else:
            # One has interval metadata while the other does not.
            return False

    @abc.abstractmethod
    def __ne__(self, other):
        raise NotImplementedError

    def _ne_(self, other):
        return not (self == other)

    @abc.abstractmethod
    def __copy__(self):
        raise NotImplementedError

    def _copy_(self):
        if self.has_interval_metadata():
            return copy.copy(self.interval_metadata)
        else:
            return None

    @abc.abstractmethod
    def __deepcopy__(self, memo):
        raise NotImplementedError

    def _deepcopy_(self, memo):
        if self.has_interval_metadata():
            return copy.deepcopy(self.interval_metadata, memo)
        else:
            return None
