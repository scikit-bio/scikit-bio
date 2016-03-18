# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import abc
import copy

import numpy as np
import pandas as pd

from skbio.util._decorator import stable


class MetadataMixin(metaclass=abc.ABCMeta):
    @property
    @stable(as_of="0.4.0")
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

        >>> from pprint import pprint
        >>> from skbio import Sequence
        >>> seq = Sequence('ACGT', metadata={'id': 'seq-id',
        ...                                  'description': 'seq description'})

        Retrieve metadata:

        >>> pprint(seq.metadata) # using pprint to display dict in sorted order
        {'description': 'seq description', 'id': 'seq-id'}

        Update metadata:

        >>> seq.metadata['id'] = 'new-id'
        >>> seq.metadata['pubmed'] = 12345
        >>> pprint(seq.metadata)
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
            raise TypeError("metadata must be a dict")
        # Shallow copy.
        self._metadata = metadata.copy()

    @metadata.deleter
    def metadata(self):
        self._metadata = None

    @abc.abstractmethod
    def __init__(self, metadata=None):
        pass

    def _init_(self, metadata=None):
        if metadata is None:
            self._metadata = None
        else:
            self.metadata = metadata

    @abc.abstractmethod
    def __eq__(self, other):
        pass

    def _eq_(self, other):
        # We're not simply comparing self.metadata to other.metadata in order
        # to avoid creating "empty" metadata representations on the objects if
        # they don't have metadata.
        if self.has_metadata() and other.has_metadata():
            if self.metadata != other.metadata:
                return False
        elif not (self.has_metadata() or other.has_metadata()):
            # Both don't have metadata.
            pass
        else:
            # One has metadata while the other does not.
            return False

        return True

    @abc.abstractmethod
    def __ne__(self, other):
        pass

    def _ne_(self, other):
        return not (self == other)

    @abc.abstractmethod
    def __copy__(self):
        pass

    def _copy_(self):
        if self.has_metadata():
            return self.metadata.copy()
        else:
            return None

    @abc.abstractmethod
    def __deepcopy__(self, memo):
        pass

    def _deepcopy_(self, memo):
        if self.has_metadata():
            return copy.deepcopy(self.metadata, memo)
        else:
            return None

    @stable(as_of="0.4.0")
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
        pass

    @property
    @stable(as_of="0.4.0")
    def positional_metadata(self):
        """``pd.DataFrame`` containing metadata along an axis.

        Notes
        -----
        This property can be set and deleted. When setting new positional
        metadata a shallow copy is made.

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
        ...     positional_metadata={'quality': [3, 3, 20, 11],
        ...                          'exons': [True, True, False, True]})
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
            has non-degenerates: True
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
        >>> seq.positional_metadata
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
                index=np.arange(self._positional_metadata_axis_len_()))
        return self._positional_metadata

    @positional_metadata.setter
    def positional_metadata(self, positional_metadata):
        try:
            # Pass copy=True to copy underlying data buffer.
            positional_metadata = pd.DataFrame(positional_metadata, copy=True)
        except pd.core.common.PandasError as e:
            raise TypeError(
                "Invalid positional metadata. Must be consumable by "
                "`pd.DataFrame` constructor. Original pandas error message: "
                "\"%s\"" % e)

        num_rows = len(positional_metadata.index)
        axis_len = self._positional_metadata_axis_len_()
        if num_rows != axis_len:
            raise ValueError(
                "Number of positional metadata values (%d) must match the "
                "positional metadata axis length (%d)."
                % (num_rows, axis_len))

        positional_metadata.reset_index(drop=True, inplace=True)
        self._positional_metadata = positional_metadata

    @positional_metadata.deleter
    def positional_metadata(self):
        self._positional_metadata = None

    @abc.abstractmethod
    def __init__(self, positional_metadata=None):
        pass

    def _init_(self, positional_metadata=None):
        if positional_metadata is None:
            self._positional_metadata = None
        else:
            self.positional_metadata = positional_metadata

    @abc.abstractmethod
    def __eq__(self, other):
        pass

    def _eq_(self, other):
        # We're not simply comparing self.positional_metadata to
        # other.positional_metadata in order to avoid creating "empty"
        # positional metadata representations on the objects if they don't have
        # positional metadata.
        if self.has_positional_metadata() and other.has_positional_metadata():
            if not self.positional_metadata.equals(other.positional_metadata):
                return False
        elif not (self.has_positional_metadata() or
                  other.has_positional_metadata()):
            # Both don't have positional metadata.
            pass
        else:
            # One has positional metadata while the other does not.
            return False

        return True

    @abc.abstractmethod
    def __ne__(self, other):
        pass

    def _ne_(self, other):
        return not (self == other)

    @abc.abstractmethod
    def __copy__(self):
        pass

    def _copy_(self):
        if self.has_positional_metadata():
            # deep=True makes a shallow copy of the underlying data buffer.
            return self.positional_metadata.copy(deep=True)
        else:
            return None

    @abc.abstractmethod
    def __deepcopy__(self, memo):
        pass

    def _deepcopy_(self, memo):
        if self.has_positional_metadata():
            return copy.deepcopy(self.positional_metadata, memo)
        else:
            return None

    @stable(as_of="0.4.0")
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
        return (self._positional_metadata is not None and
                len(self.positional_metadata.columns) > 0)
