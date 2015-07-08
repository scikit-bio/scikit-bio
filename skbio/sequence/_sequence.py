# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range
from future.utils import viewitems
import six

import itertools
import math
import re
import collections
import copy
import numbers
import textwrap
from contextlib import contextmanager

import numpy as np
from scipy.spatial.distance import hamming

import pandas as pd

from skbio._base import SkbioObject
from skbio.sequence._base import ElasticLines
from skbio.util._misc import chunk_str


class Sequence(collections.Sequence, SkbioObject):
    """Store biological sequence data and optional associated metadata.

    ``Sequence`` objects do not enforce an alphabet and are thus the most
    generic objects for storing biological sequence data. Subclasses ``DNA``,
    ``RNA``, and ``Protein`` enforce the IUPAC character set [1]_ for, and
    provide operations specific to, each respective molecule type.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
        Characters representing the biological sequence itself.
    metadata : dict, optional
        Arbitrary metadata which applies to the entire sequence. A shallow copy
        of the ``dict`` will be made (see Examples section below for details).
    positional_metadata : pd.DataFrame consumable, optional
        Arbitrary per-character metadata (e.g., sequence read quality
        scores). Must be able to be passed directly to ``pd.DataFrame``
        constructor. Each column of metadata must be the same length as the
        biological sequence. A shallow copy of the positional metadata will be
        made if necessary (see Examples section below for details).

    Attributes
    ----------
    values
    metadata
    positional_metadata

    See Also
    --------
    DNA
    RNA
    Protein

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    Examples
    --------
    >>> from pprint import pprint
    >>> from skbio import Sequence

    **Creating sequences:**

    Create a sequence without any metadata:

    >>> seq = Sequence('GGUCGUGAAGGA')
    >>> seq
    Sequence
    ---------------
    Stats:
        length: 12
    ---------------
    0 GGUCGUGAAG GA

    Create a sequence with metadata and positional metadata:

    >>> metadata = {'id':'seq-id', 'desc':'seq desc', 'authors': ['Alice']}
    >>> positional_metadata = {'quality': [3, 3, 4, 10],
    ...                        'exons': [True, True, False, True]}
    >>> seq = Sequence('ACGT', metadata=metadata,
    ...                positional_metadata=positional_metadata)
    >>> seq
    Sequence
    -----------------------------
    Metadata:
        'authors': <type 'list'>
        'desc': 'seq desc'
        'id': 'seq-id'
    Positional metadata:
        'exons': <dtype: bool>
        'quality': <dtype: int64>
    Stats:
        length: 4
    -----------------------------
    0 ACGT

    **Retrieving sequence metadata:**

    Retrieve metadata:

    >>> pprint(seq.metadata) # using pprint to display dict in sorted order
    {'authors': ['Alice'], 'desc': 'seq desc', 'id': 'seq-id'}

    Retrieve positional metadata:

    >>> seq.positional_metadata
       exons  quality
    0   True        3
    1   True        3
    2  False        4
    3   True       10

    **Updating sequence metadata:**

    .. warning:: Be aware that a shallow copy of ``metadata`` and
       ``positional_metadata`` is made for performance. Since a deep copy is
       not made, changes made to mutable Python objects stored as metadata may
       affect the metadata of other ``Sequence`` objects or anything else that
       shares a reference to the object. The following examples illustrate this
       behavior.

    First, let's create a sequence and update its metadata:

    >>> metadata = {'id':'seq-id', 'desc':'seq desc', 'authors': ['Alice']}
    >>> seq = Sequence('ACGT', metadata=metadata)
    >>> seq.metadata['id'] = 'new-id'
    >>> seq.metadata['pubmed'] = 12345
    >>> pprint(seq.metadata)
    {'authors': ['Alice'], 'desc': 'seq desc', 'id': 'new-id', 'pubmed': 12345}

    Note that the original metadata dictionary (stored in variable
    ``metadata``) hasn't changed because a shallow copy was made:

    >>> pprint(metadata)
    {'authors': ['Alice'], 'desc': 'seq desc', 'id': 'seq-id'}
    >>> seq.metadata == metadata
    False

    Note however that since only a *shallow* copy was made, updates to mutable
    objects will also change the original metadata dictionary:

    >>> seq.metadata['authors'].append('Bob')
    >>> seq.metadata['authors']
    ['Alice', 'Bob']
    >>> metadata['authors']
    ['Alice', 'Bob']

    This behavior can also occur when manipulating a sequence that has been
    derived from another sequence:

    >>> subseq = seq[1:3]
    >>> subseq
    Sequence
    ----------------------------
    Metadata:
        'authors': <type 'list'>
        'desc': 'seq desc'
        'id': 'new-id'
        'pubmed': 12345
    Stats:
        length: 2
    ----------------------------
    0 CG
    >>> pprint(subseq.metadata)
    {'authors': ['Alice', 'Bob'],
     'desc': 'seq desc',
     'id': 'new-id',
     'pubmed': 12345}

    The subsequence has inherited the metadata of its parent sequence. If we
    update the subsequence's author list, we see the changes propagated in the
    parent sequence and original metadata dictionary:

    >>> subseq.metadata['authors'].append('Carol')
    >>> subseq.metadata['authors']
    ['Alice', 'Bob', 'Carol']
    >>> seq.metadata['authors']
    ['Alice', 'Bob', 'Carol']
    >>> metadata['authors']
    ['Alice', 'Bob', 'Carol']

    The behavior for updating positional metadata is similar. Let's create a
    new sequence with positional metadata that is already stored in a
    ``pd.DataFrame``:

    >>> positional_metadata = pd.DataFrame(
    ...     {'quality': [3, 3, 4, 10], 'list': [[], [], [], []]})
    >>> seq = Sequence('ACGT', positional_metadata=positional_metadata)
    >>> seq
    Sequence
    -----------------------------
    Positional metadata:
        'list': <dtype: object>
        'quality': <dtype: int64>
    Stats:
        length: 4
    -----------------------------
    0 ACGT
    >>> seq.positional_metadata
      list  quality
    0   []        3
    1   []        3
    2   []        4
    3   []       10

    Now let's update the sequence's positional metadata by adding a new column
    and changing a value in another column:

    >>> seq.positional_metadata['gaps'] = [False, False, False, False]
    >>> seq.positional_metadata.loc[0, 'quality'] = 999
    >>> seq.positional_metadata
      list  quality   gaps
    0   []      999  False
    1   []        3  False
    2   []        4  False
    3   []       10  False

    Note that the original positional metadata (stored in variable
    ``positional_metadata``) hasn't changed because a shallow copy was made:

    >>> positional_metadata
      list  quality
    0   []        3
    1   []        3
    2   []        4
    3   []       10
    >>> seq.positional_metadata.equals(positional_metadata)
    False

    Next let's create a sequence that has been derived from another sequence:

    >>> subseq = seq[1:3]
    >>> subseq
    Sequence
    -----------------------------
    Positional metadata:
        'list': <dtype: object>
        'quality': <dtype: int64>
        'gaps': <dtype: bool>
    Stats:
        length: 2
    -----------------------------
    0 CG
    >>> subseq.positional_metadata
      list  quality   gaps
    0   []        3  False
    1   []        4  False

    As described above for metadata, since only a *shallow* copy was made of
    the positional metadata, updates to mutable objects will also change the
    parent sequence's positional metadata and the original positional metadata
    ``pd.DataFrame``:

    >>> subseq.positional_metadata.loc[0, 'list'].append('item')
    >>> subseq.positional_metadata
         list  quality   gaps
    0  [item]        3  False
    1      []        4  False
    >>> seq.positional_metadata
         list  quality   gaps
    0      []      999  False
    1  [item]        3  False
    2      []        4  False
    3      []       10  False
    >>> positional_metadata
         list  quality
    0      []        3
    1  [item]        3
    2      []        4
    3      []       10

    """
    default_write_format = 'fasta'
    __hash__ = None

    @property
    def values(self):
        """Array containing underlying sequence characters.

        Notes
        -----
        This property is not writeable.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('AACGA')
        >>> s.values # doctest: +NORMALIZE_WHITESPACE
        array(['A', 'A', 'C', 'G', 'A'],
              dtype='|S1')

        """
        return self._bytes.view('|S1')

    @property
    def metadata(self):
        """``dict`` containing metadata which applies to the entire sequence.

        Notes
        -----
        This property can be set and deleted.

        Examples
        --------
        >>> from pprint import pprint
        >>> from skbio import Sequence

        Create a sequence with metadata:

        >>> s = Sequence('ACGTACGTACGTACGT',
        ...              metadata={'id': 'seq-id',
        ...                        'description': 'seq description'})
        >>> s
        Sequence
        ------------------------------------
        Metadata:
            'description': 'seq description'
            'id': 'seq-id'
        Stats:
            length: 16
        ------------------------------------
        0 ACGTACGTAC GTACGT

        Retrieve metadata:

        >>> pprint(s.metadata) # using pprint to display dict in sorted order
        {'description': 'seq description', 'id': 'seq-id'}

        Update metadata:

        >>> s.metadata['id'] = 'new-id'
        >>> s.metadata['pubmed'] = 12345
        >>> pprint(s.metadata)
        {'description': 'seq description', 'id': 'new-id', 'pubmed': 12345}

        Set metadata:

        >>> s.metadata = {'abc': 123}
        >>> s.metadata
        {'abc': 123}

        Delete metadata:

        >>> s.has_metadata()
        True
        >>> del s.metadata
        >>> s.metadata
        {}
        >>> s.has_metadata()
        False

        """
        if self._metadata is None:
            # not using setter to avoid copy
            self._metadata = {}
        return self._metadata

    @metadata.setter
    def metadata(self, metadata):
        if not isinstance(metadata, dict):
            raise TypeError("metadata must be a dict")
        # shallow copy
        self._metadata = metadata.copy()

    @metadata.deleter
    def metadata(self):
        self._metadata = None

    @property
    def positional_metadata(self):
        """``pd.DataFrame`` containing metadata on a per-character basis.

        Notes
        -----
        This property can be set and deleted.

        Examples
        --------
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
            # not using setter to avoid copy
            self._positional_metadata = pd.DataFrame(
                index=np.arange(len(self)))
        return self._positional_metadata

    @positional_metadata.setter
    def positional_metadata(self, positional_metadata):
        try:
            # copy=True to copy underlying data buffer
            positional_metadata = pd.DataFrame(positional_metadata, copy=True)
        except pd.core.common.PandasError as e:
            raise TypeError('Positional metadata invalid. Must be consumable '
                            'by pd.DataFrame. Original pandas error message: '
                            '"%s"' % e)

        num_rows = len(positional_metadata.index)
        if num_rows != len(self):
            raise ValueError(
                "Number of positional metadata values (%d) must match the "
                "number of characters in the sequence (%d)." %
                (num_rows, len(self)))

        positional_metadata.reset_index(drop=True, inplace=True)
        self._positional_metadata = positional_metadata

    @positional_metadata.deleter
    def positional_metadata(self):
        self._positional_metadata = None

    @property
    def _string(self):
        return self._bytes.tostring()

    def __init__(self, sequence, metadata=None,
                 positional_metadata=None):

        if isinstance(sequence, np.ndarray):
            if sequence.dtype == np.uint8:
                self._set_bytes_contiguous(sequence)
            elif sequence.dtype == '|S1':
                sequence = sequence.view(np.uint8)
                # Guarantee the sequence is an array (might be scalar before
                # this).
                if sequence.shape == ():
                    sequence = np.array([sequence], dtype=np.uint8)
                self._set_bytes_contiguous(sequence)
            else:
                raise TypeError(
                    "Can only create sequence from numpy.ndarray of dtype "
                    "np.uint8 or '|S1'. Invalid dtype: %s" %
                    sequence.dtype)
        elif isinstance(sequence, Sequence):
            # we're not simply accessing sequence.metadata in order to avoid
            # creating "empty" metadata representations on both sequence
            # objects if they don't have metadata. same strategy is used below
            # for positional metadata
            if metadata is None and sequence.has_metadata():
                metadata = sequence.metadata
            if (positional_metadata is None and
                    sequence.has_positional_metadata()):
                positional_metadata = sequence.positional_metadata
            sequence = sequence._bytes

            self._owns_bytes = False

            self._set_bytes(sequence)

        else:
            # Python 3 will not raise a UnicodeEncodeError so we force it by
            # encoding it as ascii
            if isinstance(sequence, six.text_type):
                sequence = sequence.encode("ascii")
            s = np.fromstring(sequence, dtype=np.uint8)

            # There are two possibilities (to our knowledge) at this point:
            # Either the sequence we were given was something string-like,
            # (else it would not have made it past fromstring), or it was a
            # numpy scalar, and so our length must be 1.
            if isinstance(sequence, np.generic) and len(s) != 1:
                raise TypeError("Can cannot create a sequence with %r" %
                                type(sequence).__name__)

            sequence = s
            self._owns_bytes = True

            self._set_bytes(sequence)

        if metadata is None:
            self._metadata = None
        else:
            self.metadata = metadata

        if positional_metadata is None:
            self._positional_metadata = None
        else:
            self.positional_metadata = positional_metadata

    def _set_bytes_contiguous(self, sequence):
        """Munge the sequence data into a numpy array of dtype uint8."""
        if not sequence.flags['C_CONTIGUOUS']:
            # numpy doesn't support views of non-contiguous arrays. Since we're
            # making heavy use of views internally, and users may also supply
            # us with a view, make sure we *always* store a contiguous array to
            # avoid hard-to-track bugs. See
            # https://github.com/numpy/numpy/issues/5716
            sequence = np.ascontiguousarray(sequence)
            self._owns_bytes = True
        else:
            self._owns_bytes = False
        self._set_bytes(sequence)

    def _set_bytes(self, sequence):
        sequence.flags.writeable = False
        self._bytes = sequence

    def __contains__(self, subsequence):
        """Determine if a subsequence is contained in the biological sequence.

        Parameters
        ----------
        subsequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            The putative subsequence.

        Returns
        -------
        bool
            Indicates whether `subsequence` is contained in the biological
            sequence.

        Raises
        ------
        TypeError
            If `subsequence` is a ``Sequence`` object with a different type
            than the biological sequence.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> 'GGU' in s
        True
        >>> 'CCC' in s
        False

        """
        return self._munge_to_bytestring(subsequence, "in") in self._string

    def __eq__(self, other):
        """Determine if the biological sequence is equal to another.

        Biological sequences are equal if they are *exactly* the same type and
        their sequence characters, metadata, and positional metadata are the
        same.

        Parameters
        ----------
        other : Sequence
            Sequence to test for equality against.

        Returns
        -------
        bool
            Indicates whether the biological sequence is equal to `other`.

        Examples
        --------
        Define two biological sequences that have the same underlying sequence
        of characters:

        >>> from skbio import Sequence
        >>> s = Sequence('ACGT')
        >>> t = Sequence('ACGT')

        The two sequences are considered equal because they are the same type,
        their underlying sequence of characters are the same, and their
        optional metadata attributes (``metadata`` and ``positional_metadata``)
        were not provided:

        >>> s == t
        True
        >>> t == s
        True

        Define another biological sequence with a different sequence of
        characters than the previous two biological sequences:

        >>> u = Sequence('ACGA')
        >>> u == t
        False

        Define a biological sequence with the same sequence of characters as
        ``u`` but with different metadata and positional metadata:

        >>> v = Sequence('ACGA', metadata={'id': 'abc'},
        ...              positional_metadata={'quality':[1, 5, 3, 3]})

        The two sequences are not considered equal because their metadata and
        positional metadata do not match:

        >>> u == v
        False

        """
        # checks ordered from least to most expensive
        if self.__class__ != other.__class__:
            return False

        # we're not simply comparing self.metadata to other.metadata in order
        # to avoid creating "empty" metadata representations on the sequence
        # objects if they don't have metadata. same strategy is used below for
        # positional metadata
        if self.has_metadata() and other.has_metadata():
            if self.metadata != other.metadata:
                return False
        elif not (self.has_metadata() or other.has_metadata()):
            # both don't have metadata
            pass
        else:
            # one has metadata while the other does not
            return False

        if self._string != other._string:
            return False

        if self.has_positional_metadata() and other.has_positional_metadata():
            if not self.positional_metadata.equals(other.positional_metadata):
                return False
        elif not (self.has_positional_metadata() or
                  other.has_positional_metadata()):
            # both don't have positional metadata
            pass
        else:
            # one has positional metadata while the other does not
            return False

        return True

    def __ne__(self, other):
        """Determine if the biological sequence is not equal to another.

        Biological sequences are not equal if they are not *exactly* the same
        type, or their sequence characters, metadata, or positional metadata
        differ.

        Parameters
        ----------
        other : Sequence
            Sequence to test for inequality against.

        Returns
        -------
        bool
            Indicates whether the biological sequence is not equal to `other`.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('ACGT')
        >>> t = Sequence('ACGT')
        >>> s != t
        False
        >>> u = Sequence('ACGA')
        >>> u != t
        True
        >>> v = Sequence('ACGA', metadata={'id': 'v'})
        >>> u != v
        True

        """
        return not (self == other)

    def __getitem__(self, indexable):
        """Slice the biological sequence.

        Parameters
        ----------
        indexable : int, slice, iterable (int and slice), 1D array_like (bool)
            The position(s) to return from the biological sequence. If
            `indexable` is an iterable of integers, these are assumed to be
            indices in the sequence to keep. If `indexable` is a 1D
            ``array_like`` of booleans, these are assumed to be the positions
            in the sequence to keep.

        Returns
        -------
        Sequence
            New biological sequence containing the position(s) specified by
            `indexable` in the current biological sequence. If quality scores
            are present, they will be sliced in the same manner and included in
            the returned biological sequence. ID and description are also
            included.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')

        Obtain a single character from the biological sequence:

        >>> s[1]
        Sequence
        -------------
        Stats:
            length: 1
        -------------
        0 G

        Obtain a slice:

        >>> s[7:]
        Sequence
        -------------
        Stats:
            length: 5
        -------------
        0 AAGGA

        Obtain characters at the following indices:

        >>> s[[3, 4, 7, 0, 3]]
        Sequence
        -------------
        Stats:
            length: 5
        -------------
        0 CGAGC

        Obtain characters at positions evaluating to `True`:

        >>> s = Sequence('GGUCG')
        >>> index = [True, False, True, 'a' is 'a', False]
        >>> s[index]
        Sequence
        -------------
        Stats:
            length: 3
        -------------
        0 GUC

        """
        if (not isinstance(indexable, np.ndarray) and
            ((not isinstance(indexable, six.string_types)) and
             hasattr(indexable, '__iter__'))):
            indexable_ = indexable
            indexable = np.asarray(indexable)

            if indexable.dtype == object:
                indexable = list(indexable_)  # TODO: Don't blow out memory

                if len(indexable) == 0:
                    # indexing with an empty list, so convert to ndarray and
                    # fall through to ndarray slicing below
                    indexable = np.asarray(indexable)
                else:
                    seq = np.concatenate(
                        list(_slices_from_iter(self._bytes, indexable)))
                    index = _as_slice_if_single_index(indexable)

                    positional_metadata = None
                    if self.has_positional_metadata():
                        pos_md_slices = list(_slices_from_iter(
                                             self.positional_metadata, index))
                        positional_metadata = pd.concat(pos_md_slices)

                    return self._to(sequence=seq,
                                    positional_metadata=positional_metadata)
        elif (isinstance(indexable, six.string_types) or
                isinstance(indexable, bool)):
            raise IndexError("Cannot index with %s type: %r" %
                             (type(indexable).__name__, indexable))

        if (isinstance(indexable, np.ndarray) and
            indexable.dtype == bool and
                len(indexable) != len(self)):
            raise IndexError("An boolean vector index must be the same length"
                             " as the sequence (%d, not %d)." %
                             (len(self), len(indexable)))

        if isinstance(indexable, np.ndarray) and indexable.size == 0:
            # convert an empty ndarray to a supported dtype for slicing a numpy
            # array
            indexable = indexable.astype(int)

        seq = self._bytes[indexable]
        positional_metadata = self._slice_positional_metadata(indexable)

        return self._to(sequence=seq, positional_metadata=positional_metadata)

    def _slice_positional_metadata(self, indexable):
        if self.has_positional_metadata():
            if _is_single_index(indexable):
                index = _single_index_to_slice(indexable)
            else:
                index = indexable
            return self.positional_metadata.iloc[index]
        else:
            return None

    def __len__(self):
        """Return the number of characters in the biological sequence.

        Returns
        -------
        int
            The length of the biological sequence.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> len(s)
        4

        """
        return self._bytes.size

    def __nonzero__(self):
        """Returns truth value (truthiness) of sequence.

        Returns
        -------
        bool
            True if length of sequence is greater than 0, else False.

        Examples
        --------
        >>> from skbio import Sequence
        >>> bool(Sequence(''))
        False
        >>> bool(Sequence('ACGT'))
        True

        """
        return len(self) > 0

    def __iter__(self):
        """Iterate over positions in the biological sequence.

        Yields
        ------
        Sequence
            Single character subsequence, one for each position in the
            sequence.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> for c in s:
        ...     str(c)
        'G'
        'G'
        'U'
        'C'

        """
        for i in range(len(self)):
            yield self[i]

    def __reversed__(self):
        """Iterate over positions in the biological sequence in reverse order.

        Yields
        ------
        Sequence
            Single character subsequence, one for each position in the
            sequence.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> for c in reversed(s):
        ...     str(c)
        'C'
        'U'
        'G'
        'G'

        """
        return iter(self[::-1])

    def __str__(self):
        """Return biological sequence characters as a string.

        Returns
        -------
        str
            Sequence characters as a string. No metadata or positional
            metadata will be included.

        See Also
        --------
        sequence

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUAAAGGA', metadata={'id':'hello'})
        >>> str(s)
        'GGUCGUAAAGGA'

        """
        return str(self._string.decode("ascii"))

    def __repr__(self):
        r"""Return a string representation of the biological sequence object.

        Representation includes:

        * sequence type
        * metadata keys and values: will display key/value if it is an
          understood type, otherwise just the type will be displayed. If it is
          an understood type whose representation is too long, just the type
          will be displayed
        * positional metadata: column names and column dtypes will be displayed
          in the order they appear in the positional metadata ``pd.DataFrame``.
          Column names (i.e., keys) follow the same display rules as metadata
          keys
        * sequence stats (e.g., length)
        * up to five lines of chunked sequence data. Each line of chunked
          sequence data displays the current position in the sequence

        Returns
        -------
        str
            String representation of the biological sequence object.

        Notes
        -----
        Subclasses can override Sequence._repr_stats to provide custom
        statistics.

        Examples
        --------
        Short sequence without metadata:

        >>> from skbio import Sequence
        >>> Sequence('ACGTAATGGATACGTAATGCA')
        Sequence
        -------------------------
        Stats:
            length: 21
        -------------------------
        0 ACGTAATGGA TACGTAATGC A

        Longer sequence displays first two lines and last two lines:

        >>> Sequence('ACGT' * 100)
        Sequence
        ---------------------------------------------------------------------
        Stats:
            length: 400
        ---------------------------------------------------------------------
        0   ACGTACGTAC GTACGTACGT ACGTACGTAC GTACGTACGT ACGTACGTAC GTACGTACGT
        60  ACGTACGTAC GTACGTACGT ACGTACGTAC GTACGTACGT ACGTACGTAC GTACGTACGT
        ...
        300 ACGTACGTAC GTACGTACGT ACGTACGTAC GTACGTACGT ACGTACGTAC GTACGTACGT
        360 ACGTACGTAC GTACGTACGT ACGTACGTAC GTACGTACGT

        Sequence with metadata and positional metadata:

        >>> metadata = {
        ...     'id': 'seq-id',
        ...     'description': 'description of the sequence, wrapping across '
        ...     'lines if it\'s too long',
        ...     'authors': ['Alice', 'Bob', 'Carol'],
        ...     'year': 2015,
        ...     'published': True
        ... }
        >>> positional_metadata = {
        ...     'quality': [3, 10, 11, 10],
        ...     'exons': [True, True, False, True]
        ... }
        >>> Sequence('ACGT', metadata=metadata,
        ...          positional_metadata=positional_metadata)
        Sequence
        ----------------------------------------------------------------------
        Metadata:
            'authors': <type 'list'>
            'description': "description of the sequence, wrapping across lines
                            if it's too long"
            'id': 'seq-id'
            'published': True
            'year': 2015
        Positional metadata:
            'exons': <dtype: bool>
            'quality': <dtype: int64>
        Stats:
            length: 4
        ----------------------------------------------------------------------
        0 ACGT

        """
        return _SequenceReprBuilder(
            seq=self,
            width=71,  # 79 for pep8, 8 space indent for docstrings
            indent=4,
            chunk_size=10).build()

    def _repr_stats(self):
        """Define statistics to display in the sequence's repr.

        Subclasses can override this method to provide type-specific
        statistics.

        This method computes a single statistic: length.

        Returns
        -------
        list
            List of tuples where each tuple represents a statistic. Each tuple
            contains exactly two ``str`` elements: the statistic's name/label,
            and the str-formatted value of the statistic. Ordering of
            statistics (i.e., list order) determines display order in the
            sequence repr.

        """
        return [('length', '%d' % len(self))]

    def __copy__(self):
        """Return a shallow copy of the biological sequence.

        See Also
        --------
        copy

        Notes
        -----
        This method is equivalent to ``seq.copy(deep=False)``.

        """
        return self.copy(deep=False)

    def __deepcopy__(self, memo):
        """Return a deep copy of the biological sequence.

        See Also
        --------
        copy

        Notes
        -----
        This method is equivalent to ``seq.copy(deep=True)``.

        """
        return self._copy(True, memo)

    def has_metadata(self):
        """Determine if the sequence contains metadata.

        Returns
        -------
        bool
            Indicates whether the sequence has metadata

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACACGACGTT')
        >>> s.has_metadata()
        False
        >>> t = DNA('ACACGACGTT', metadata={'id': 'seq-id'})
        >>> t.has_metadata()
        True

        """
        return self._metadata is not None and bool(self.metadata)

    def has_positional_metadata(self):
        """Determine if the sequence contains positional metadata.

        Returns
        -------
        bool
            Indicates whether the sequence has positional metadata

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACACGACGTT')
        >>> s.has_positional_metadata()
        False
        >>> t = DNA('ACACGACGTT', positional_metadata={'quality': range(10)})
        >>> t.has_positional_metadata()
        True

        """
        return (self._positional_metadata is not None and
                len(self.positional_metadata.columns) > 0)

    def copy(self, deep=False):
        """Return a copy of the biological sequence.

        Parameters
        ----------
        deep : bool, optional
            Perform a deep copy. If ``False``, perform a shallow copy.

        Returns
        -------
        Sequence
            Copy of the biological sequence.

        Notes
        -----
        Since sequence objects can share the same underlying immutable sequence
        data (or pieces of it), this method can be used to create a sequence
        object with its own copy of the sequence data so that the original
        sequence data can be garbage-collected.

        Examples
        --------
        Create a sequence:

        >>> from pprint import pprint
        >>> from skbio import Sequence
        >>> seq = Sequence('ACGT',
        ...                metadata={'id': 'seq-id', 'authors': ['Alice']},
        ...                positional_metadata={'quality': [7, 10, 8, 5],
        ...                                     'list': [[], [], [], []]})

        Make a shallow copy of the sequence:

        >>> seq_copy = seq.copy()
        >>> seq_copy == seq
        True

        Setting new references in the copied sequence's metadata doesn't affect
        the original sequence's metadata:

        >>> seq_copy.metadata['id'] = 'new-id'
        >>> pprint(seq_copy.metadata)
        {'authors': ['Alice'], 'id': 'new-id'}
        >>> pprint(seq.metadata)
        {'authors': ['Alice'], 'id': 'seq-id'}

        The same applies to the sequence's positional metadata:

        >>> seq_copy.positional_metadata.loc[0, 'quality'] = 999
        >>> seq_copy.positional_metadata
          list  quality
        0   []      999
        1   []       10
        2   []        8
        3   []        5
        >>> seq.positional_metadata
          list  quality
        0   []        7
        1   []       10
        2   []        8
        3   []        5

        Since only a *shallow* copy was made, updates to mutable objects stored
        as metadata affect the original sequence's metadata:

        >>> seq_copy.metadata['authors'].append('Bob')
        >>> pprint(seq_copy.metadata)
        {'authors': ['Alice', 'Bob'], 'id': 'new-id'}
        >>> pprint(seq.metadata)
        {'authors': ['Alice', 'Bob'], 'id': 'seq-id'}

        The same applies to the sequence's positional metadata:

        >>> seq_copy.positional_metadata.loc[0, 'list'].append(1)
        >>> seq_copy.positional_metadata
          list  quality
        0  [1]      999
        1   []       10
        2   []        8
        3   []        5
        >>> seq.positional_metadata
          list  quality
        0  [1]        7
        1   []       10
        2   []        8
        3   []        5

        Perform a deep copy to avoid this behavior:

        >>> seq_deep_copy = seq.copy(deep=True)

        Updates to mutable objects no longer affect the original sequence's
        metadata:

        >>> seq_deep_copy.metadata['authors'].append('Carol')
        >>> pprint(seq_deep_copy.metadata)
        {'authors': ['Alice', 'Bob', 'Carol'], 'id': 'seq-id'}
        >>> pprint(seq.metadata)
        {'authors': ['Alice', 'Bob'], 'id': 'seq-id'}

        Nor its positional metadata:

        >>> seq_deep_copy.positional_metadata.loc[0, 'list'].append(2)
        >>> seq_deep_copy.positional_metadata
             list  quality
        0  [1, 2]        7
        1      []       10
        2      []        8
        3      []        5
        >>> seq.positional_metadata
          list  quality
        0  [1]        7
        1   []       10
        2   []        8
        3   []        5

        """
        return self._copy(deep, {})

    def _copy(self, deep, memo):
        # strategy: copy the sequence without metadata first, then set metadata
        # attributes with copies. we take this approach instead of simply
        # passing the metadata through the Sequence constructor because we
        # don't want to copy twice (this could happen when deep=True, where we
        # deep copy here and then shallow copy in the Sequence constructor). we
        # also directly set the private metadata attributes instead of using
        # their public setters to avoid an unnecessary copy

        # we don't make a distinction between deep vs. shallow copy of bytes
        # because dtype=np.uint8. we only need to make the distinction when
        # dealing with object dtype
        bytes = np.copy(self._bytes)

        seq_copy = self._constructor(sequence=bytes, metadata=None,
                                     positional_metadata=None)

        if self.has_metadata():
            metadata = self.metadata
            if deep:
                metadata = copy.deepcopy(metadata, memo)
            else:
                metadata = metadata.copy()
            seq_copy._metadata = metadata

        if self.has_positional_metadata():
            positional_metadata = self.positional_metadata
            if deep:
                positional_metadata = copy.deepcopy(positional_metadata, memo)
            else:
                # deep=True makes a shallow copy of the underlying data buffer
                positional_metadata = positional_metadata.copy(deep=True)
            seq_copy._positional_metadata = positional_metadata

        return seq_copy

    def count(self, subsequence, start=None, end=None):
        """Count occurrences of a subsequence in the biological sequence.

        Parameters
        ----------
        subsequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Subsequence to count occurrences of.
        start : int, optional
            The position at which to start counting (inclusive).
        end : int, optional
            The position at which to stop counting (exclusive).

        Returns
        -------
        int
            Number of occurrences of `subsequence` in the biological sequence.

        Raises
        ------
        ValueError
            If `subsequence` is of length 0.
        TypeError
            If `subsequence` is a ``Sequence`` object with a different type
            than the biological sequence.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUCG')
        >>> s.count('G')
        3
        >>> s.count('GG')
        1
        >>> s.count('T')
        0
        >>> s.count('G', 2, 5)
        1

        """
        if len(subsequence) == 0:
            raise ValueError("`count` is not defined for empty subsequences.")

        return self._string.count(
            self._munge_to_bytestring(subsequence, "count"), start, end)

    def index(self, subsequence, start=None, end=None):
        """Find position where subsequence first occurs in the sequence.

        Parameters
        ----------
        subsequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Subsequence to search for in the biological sequence.
        start : int, optional
            The position at which to start searching (inclusive).
        end : int, optional
            The position at which to stop searching (exclusive).

        Returns
        -------
        int
            Position where `subsequence` first occurs in the biological
            sequence.

        Raises
        ------
        ValueError
            If `subsequence` is not present in the biological sequence.
        TypeError
            If `subsequence` is a ``Sequence`` object with a different type
            than the biological sequence.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('ACACGACGTT-')
        >>> s.index('ACG')
        2

        """
        try:
            return self._string.index(
                self._munge_to_bytestring(subsequence, "index"), start, end)
        except ValueError:
            raise ValueError(
                "%r is not present in %r." % (subsequence, self))

    def distance(self, other, metric=None):
        """Compute the distance to another sequence.

        Parameters
        ----------
        other : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Sequence to compute the distance to.
        metric : function, optional
            Function used to compute the distance between the biological
            sequence and `other`. If ``None`` (the default),
            ``scipy.spatial.distance.hamming`` will be used.

        Returns
        -------
        float
            Distance between the biological sequence and `other`.

        Raises
        ------
        ValueError
            If the sequences are not the same length when `metric` is ``None``
            (i.e., `metric` is ``scipy.spatial.distance.hamming``). This is
            only checked when using this metric, as equal length is not a
            requirement of all sequence distance metrics. In general, the
            metric itself should test and give an informative error message,
            but the message from ``scipy.spatial.distance.hamming`` is somewhat
            cryptic (as of this writing), and it's the default metric, so we
            explicitly do this check here. This metric-specific check will be
            removed from this method when the ``skbio.sequence.stats`` module
            is created (track progress on issue #913).
        TypeError
            If `other` is a ``Sequence`` object with a different type than the
            biological sequence.

        See Also
        --------
        fraction_diff
        fraction_same
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('AGUC')
        >>> s.distance(t)
        0.25
        >>> def custom_dist(s1, s2): return 0.42
        >>> s.distance(t, custom_dist)
        0.42

        """
        # TODO refactor this method to accept a name (string) of the distance
        # metric to apply and accept **kwargs
        other = self._munge_to_sequence(other, 'distance')
        if metric is None:
            # Hamming requires equal length sequences. We are checking this
            # here because the error you would get otherwise is cryptic.
            if len(self) != len(other):
                raise ValueError(
                    "Sequences do not have equal length. "
                    "Hamming distances can only be computed between "
                    "sequences of equal length.")
            metric = hamming
        return float(metric(self.values, other.values))

    def matches(self, other):
        """Find positions that match with another sequence.

        Parameters
        ----------
        other : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Sequence to compare to.

        Returns
        -------
        1D np.ndarray (bool)
            Boolean vector where ``True`` at position ``i`` indicates a match
            between the sequences at their positions ``i``.

        Raises
        ------
        ValueError
            If the sequences are not the same length.
        TypeError
            If `other` is a ``Sequence`` object with a different type than the
            biological sequence.

        See Also
        --------
        mismatches

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('GAUU')
        >>> s.matches(t)
        array([ True, False,  True, False], dtype=bool)

        """
        other = self._munge_to_sequence(other, 'matches/mismatches')
        if len(self) != len(other):
            raise ValueError("Match and mismatch vectors can only be "
                             "generated from equal length sequences.")
        return self._bytes == other._bytes

    def mismatches(self, other):
        """Find positions that do not match with another sequence.

        Parameters
        ----------
        other : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Sequence to compare to.

        Returns
        -------
        1D np.ndarray (bool)
            Boolean vector where ``True`` at position ``i`` indicates a
            mismatch between the sequences at their positions ``i``.

        Raises
        ------
        ValueError
            If the sequences are not the same length.
        TypeError
            If `other` is a ``Sequence`` object with a different type than the
            biological sequence.

        See Also
        --------
        matches

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('GAUU')
        >>> s.mismatches(t)
        array([False,  True, False,  True], dtype=bool)

        """
        return np.invert(self.matches(other))

    def match_frequency(self, other, relative=False):
        """Return count of positions that are the same between two sequences.

        Parameters
        ----------
        other : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Sequence to compare to.
        relative : bool, optional
            If ``True``, return the relative frequency of matches instead of
            the count.

        Returns
        -------
        int or float
            Number of positions that are the same between the sequences. This
            will be an ``int`` if `relative` is ``False`` and a ``float``
            if `relative` is ``True``.

        Raises
        ------
        ValueError
            If the sequences are not the same length.
        TypeError
            If `other` is a ``Sequence`` object with a different type than the
            biological sequence.

        See Also
        --------
        mismatch_frequency
        matches
        mismatches
        distance

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('AGUC')
        >>> s.match_frequency(t)
        3
        >>> s.match_frequency(t, relative=True)
        0.75

        """
        if relative:
            return float(self.matches(other).mean())
        else:
            return int(self.matches(other).sum())

    def mismatch_frequency(self, other, relative=False):
        """Return count of positions that differ between two sequences.

        Parameters
        ----------
        other : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Sequence to compare to.
        relative : bool, optional
            If ``True``, return the relative frequency of mismatches instead of
            the count.

        Returns
        -------
        int or float
            Number of positions that differ between the sequences. This will be
            an ``int`` if `relative` is ``False`` and a ``float``
            if `relative` is ``True``.

        Raises
        ------
        ValueError
            If the sequences are not the same length.
        TypeError
            If `other` is a ``Sequence`` object with a different type than the
            biological sequence.

        See Also
        --------
        match_frequency
        matches
        mismatches
        distance

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('AGUC')
        >>> s.mismatch_frequency(t)
        1
        >>> s.mismatch_frequency(t, relative=True)
        0.25

        """
        if relative:
            return float(self.mismatches(other).mean())
        else:
            return int(self.mismatches(other).sum())

    def iter_kmers(self, k, overlap=True):
        """Generate kmers of length `k` from the biological sequence.

        Parameters
        ----------
        k : int
            The kmer length.
        overlap : bool, optional
            Defines whether the kmers should be overlapping or not.

        Yields
        ------
        Sequence
            kmer of length `k` contained in the biological sequence.

        Raises
        ------
        ValueError
            If `k` is less than 1.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('ACACGACGTT')
        >>> for kmer in s.iter_kmers(4, overlap=False):
        ...     str(kmer)
        'ACAC'
        'GACG'
        >>> for kmer in s.iter_kmers(3, overlap=True):
        ...     str(kmer)
        'ACA'
        'CAC'
        'ACG'
        'CGA'
        'GAC'
        'ACG'
        'CGT'
        'GTT'

        """
        if k < 1:
            raise ValueError("k must be greater than 0.")

        if overlap:
            step = 1
            count = len(self) - k + 1
        else:
            step = k
            count = len(self) // k

        if self.has_positional_metadata():
            for i in range(0, len(self) - k + 1, step):
                yield self[i:i+k]
        # Optimized path when no positional metadata
        else:
            kmers = np.lib.stride_tricks.as_strided(
                self._bytes, shape=(k, count), strides=(1, step)).T
            for s in kmers:
                yield self._to(sequence=s)

    def kmer_frequencies(self, k, overlap=True, relative=False):
        """Return counts of words of length `k` from the biological sequence.

        Parameters
        ----------
        k : int
            The word length.
        overlap : bool, optional
            Defines whether the kmers should be overlapping or not.
        relative : bool, optional
            If ``True``, return the relative frequency of each kmer instead of
            its count.

        Returns
        -------
        collections.Counter or collections.defaultdict
            Frequencies of words of length `k` contained in the biological
            sequence. This will be a ``collections.Counter`` if `relative` is
            ``False`` and a ``collections.defaultdict`` if `relative` is
            ``True``.

        Raises
        ------
        ValueError
            If `k` is less than 1.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('ACACATTTATTA')
        >>> s.kmer_frequencies(3, overlap=False)
        Counter({'TTA': 2, 'ACA': 1, 'CAT': 1})
        >>> s.kmer_frequencies(3, relative=True, overlap=False)
        defaultdict(<type 'float'>, {'ACA': 0.25, 'TTA': 0.5, 'CAT': 0.25})

        """
        kmers = self.iter_kmers(k, overlap=overlap)
        freqs = collections.Counter((str(seq) for seq in kmers))

        if relative:
            if overlap:
                num_kmers = len(self) - k + 1
            else:
                num_kmers = len(self) // k

            relative_freqs = collections.defaultdict(float)
            for kmer, count in viewitems(freqs):
                relative_freqs[kmer] = count / num_kmers
            freqs = relative_freqs

        return freqs

    def find_with_regex(self, regex, ignore=None):
        """Generate slices for patterns matched by a regular expression.

        Parameters
        ----------
        regex : str or regular expression object
            String to be compiled into a regular expression, or a pre-
            compiled regular expression object (e.g., from calling
            ``re.compile``).
        ignore : 1D array_like (bool) or iterable (slices or ints), optional
            Indicate the positions to ignore when matching.

        Yields
        ------
        slice
            Location where the regular expression matched.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('AATATACCGGTTATAA')
        >>> for match in s.find_with_regex('(TATA+)'):
        ...     match
        ...     str(s[match])
        slice(2, 6, None)
        'TATA'
        slice(11, 16, None)
        'TATAA'

        """
        if isinstance(regex, six.string_types):
            regex = re.compile(regex)

        lookup = np.arange(len(self))
        if ignore is None:
            string = str(self)
        else:
            ignore = self._munge_to_index_array(ignore)
            lookup = np.delete(lookup, ignore)
            string = str(self[lookup])

        for match in regex.finditer(string):
            # We start at 1 because we don't want the group that contains all
            # other groups.
            for g in range(1, len(match.groups())+1):
                yield slice(lookup[match.start(g)],
                            lookup[match.end(g) - 1] + 1)

    def iter_contiguous(self, included, min_length=1, invert=False):
        """Yield contiguous subsequences based on `included`.

        Parameters
        ----------
        included : 1D array_like (bool) or iterable (slices or ints)
            `included` is transformed into a flat boolean vector where each
            position will either be included or skipped. All contiguous
            included positions will be yielded as a single region.
        min_length : int, optional
            The minimum length of a subsequence for it to be yielded.
            Default is 1.
        invert : bool, optional
            Whether to invert `included` such that it describes what should be
            skipped instead of included. Default is False.

        Yields
        ------
        Sequence
            Contiguous subsequence as indicated by `included`.

        Notes
        -----
        If slices provide adjacent ranges, then they will be considered the
        same contiguous subsequence.

        Examples
        --------
        Here we use `iter_contiguous` to find all of the contiguous ungapped
        sequences using a boolean vector derived from our DNA sequence.

        >>> from skbio import DNA
        >>> s = DNA('AAA--TT-CCCC-G-')
        >>> no_gaps = ~s.gaps()
        >>> for ungapped_subsequence in s.iter_contiguous(no_gaps,
        ...                                               min_length=2):
        ...     print(ungapped_subsequence)
        AAA
        TT
        CCCC

        Note how the last potential subsequence was skipped because it would
        have been smaller than our `min_length` which was set to 2.

        We can also use `iter_contiguous` on a generator of slices as is
        produced by `find_motifs` (and `find_with_regex`).

        >>> from skbio import Protein
        >>> s = Protein('ACDFNASANFTACGNPNRTESL')
        >>> for subseq in s.iter_contiguous(s.find_motifs('N-glycosylation')):
        ...     print(subseq)
        NASANFTA
        NRTE

        Note how the first subsequence contains two N-glycosylation sites. This
        happened because they were contiguous.

        """
        idx = self._munge_to_index_array(included)
        if invert:
            idx = np.delete(np.arange(len(self)), idx)

        # Adapted from http://stackoverflow.com/a/7353335/579416
        for contig in np.split(idx, np.where(np.diff(idx) != 1)[0] + 1):
            r = self[contig]
            if len(r) >= min_length:
                yield r

    def _to(self, sequence=None, metadata=None, positional_metadata=None):
        """Return a copy of the current biological sequence.

        Returns a copy of the current biological sequence, optionally with
        updated attributes specified as keyword arguments.

        Arguments are the same as those passed to the ``Sequence`` constructor.
        The returned copy will have its attributes updated based on the
        arguments. If an attribute is missing, the copy will keep the same
        attribute as the current biological sequence. Valid attribute names
        are `'sequence'`, `'metadata'`, and `'positional_metadata'`. Default
        behavior is to return a copy of the current biological sequence
        without changing any attributes.

        Parameters
        ----------
        sequence : optional
        metadata : optional
        positional_metadata : optional

        Returns
        -------
        Sequence
            Copy of the current biological sequence, optionally with updated
            attributes based on arguments. Will be the same type as the current
            biological sequence (`self`).

        Notes
        -----
        By default, `metadata` and `positional_metadata` are shallow-copied and
        the reference to `sequence` is used (without copying) for efficiency
        since `sequence` is immutable. This differs from the behavior of
        `Sequence.copy`, which will actually copy `sequence`.

        This method is the preferred way of creating new instances from an
        existing biological sequence, instead of calling
        ``self.__class__(...)``, as the latter can be error-prone (e.g.,
        it's easy to forget to propagate attributes to the new instance).

        """
        if sequence is None:
            sequence = self._bytes
        if metadata is None:
            metadata = self._metadata
        if positional_metadata is None:
            positional_metadata = self._positional_metadata
        return self._constructor(sequence=sequence, metadata=metadata,
                                 positional_metadata=positional_metadata)

    def _constructor(self, **kwargs):
        return self.__class__(**kwargs)

    def _munge_to_index_array(self, sliceable):
        """Return an index array from something isomorphic to a boolean vector.

        """
        if isinstance(sliceable, six.string_types):
            if sliceable in self.positional_metadata:
                if self.positional_metadata[sliceable].dtype == np.bool:
                    sliceable = self.positional_metadata[sliceable]
                else:
                    raise TypeError("Column '%s' in positional metadata does "
                                    "not correspond to a boolean vector" %
                                    sliceable)
            else:
                raise ValueError("No positional metadata associated with key "
                                 "'%s'" % sliceable)

        if not hasattr(sliceable, 'dtype') or (hasattr(sliceable, 'dtype') and
                                               sliceable.dtype == 'object'):
            sliceable = tuple(sliceable)
            bool_mode = False
            int_mode = False
            for s in sliceable:
                if isinstance(s, (bool, np.bool_)):
                    bool_mode = True
                elif isinstance(s, (slice, int, np.signedinteger)) or (
                        hasattr(s, 'dtype') and s.dtype != np.bool):
                    int_mode = True
                else:
                    raise TypeError("Invalid type in iterable: %s, must be one"
                                    " of {bool, int, slice, np.signedinteger}"
                                    % s.__class__.__name__)
            if bool_mode and int_mode:
                raise TypeError("Cannot provide iterable of both bool and"
                                " int.")
            sliceable = np.r_[sliceable]

        if sliceable.dtype == np.bool:
            if sliceable.size != len(self):
                raise ValueError("Boolean array (%d) does not match length of"
                                 " sequence (%d)."
                                 % (sliceable.size, len(self)))
            normalized, = np.where(sliceable)
        else:
            normalized = np.bincount(sliceable)
            if np.any(normalized > 1):
                raise ValueError("Overlapping index regions are not allowed.")

            normalized, = np.where(normalized)
            if np.any(normalized != sliceable):
                raise ValueError("Index regions are out of order.")

        return normalized

    def _munge_to_sequence(self, other, method):
        if isinstance(other, Sequence):
            if type(other) != type(self):
                raise TypeError("Cannot use %s and %s together with `%s`" %
                                (self.__class__.__name__,
                                 other.__class__.__name__, method))
            else:
                return other

        # We don't use self.__class__ or self._constructor here because we want
        # to construct the most general type of Sequence object in order to
        # avoid validation errors.
        return Sequence(other)

    def _munge_to_bytestring(self, other, method):
        if type(other) is bytes:
            return other
        elif isinstance(other, six.string_types):
            return other.encode('ascii')
        else:
            return self._munge_to_sequence(other, method)._string

    @contextmanager
    def _byte_ownership(self):
        if not self._owns_bytes:
            self._bytes = self._bytes.copy()
            self._owns_bytes = True

        self._bytes.flags.writeable = True
        yield
        self._bytes.flags.writeable = False


def _single_index_to_slice(start_index):
    end_index = None if start_index == -1 else start_index+1
    return slice(start_index, end_index)


def _is_single_index(index):
    return (isinstance(index, numbers.Integral) and
            not isinstance(index, bool))


def _as_slice_if_single_index(indexable):
    if _is_single_index(indexable):
        return _single_index_to_slice(indexable)
    else:
        return indexable


def _slices_from_iter(array, indexables):
    for i in indexables:
        if isinstance(i, slice):
            pass
        elif _is_single_index(i):
            i = _single_index_to_slice(i)
        else:
            raise IndexError("Cannot slice sequence from iterable "
                             "containing %r." % i)

        yield array[i]


class _SequenceReprBuilder(object):
    """Build a ``Sequence`` repr.

    Parameters
    ----------
    seq : Sequence
        Sequence to repr.
    width : int
        Maximum width of the repr.
    indent : int
        Number of spaces to use for indented lines.
    chunk_size: int
        Number of characters in each chunk of a sequence.

    """
    def __init__(self, seq, width, indent, chunk_size):
        self._seq = seq
        self._width = width
        self._indent = ' ' * indent
        self._chunk_size = chunk_size

    def build(self):
        lines = ElasticLines()

        cls_name = self._seq.__class__.__name__
        lines.add_line(cls_name)
        lines.add_separator()

        if self._seq.has_metadata():
            lines.add_line('Metadata:')
            # Python 3 doesn't allow sorting of mixed types so we can't just
            # use sorted() on the metadata keys. Sort first by type then sort
            # by value within each type.
            for key in self._sorted_keys_grouped_by_type(self._seq.metadata):
                value = self._seq.metadata[key]
                lines.add_lines(self._format_metadata_key_value(key, value))

        if self._seq.has_positional_metadata():
            lines.add_line('Positional metadata:')
            for key in self._seq.positional_metadata.columns.values.tolist():
                dtype = self._seq.positional_metadata[key].dtype
                lines.add_lines(
                    self._format_positional_metadata_column(key, dtype))

        lines.add_line('Stats:')
        for label, value in self._seq._repr_stats():
            lines.add_line('%s%s: %s' % (self._indent, label, value))
        lines.add_separator()

        num_lines, num_chars, column_width = self._find_optimal_seq_chunking()

        # display entire sequence if we can, else display the first two and
        # last two lines separated by ellipsis
        if num_lines <= 5:
            lines.add_lines(self._format_chunked_seq(
                range(num_lines), num_chars, column_width))
        else:
            lines.add_lines(self._format_chunked_seq(
                range(2), num_chars, column_width))
            lines.add_line('...')
            lines.add_lines(self._format_chunked_seq(
                range(num_lines - 2, num_lines), num_chars, column_width))

        return lines.to_str()

    def _sorted_keys_grouped_by_type(self, dict_):
        """Group keys within a dict by their type and sort within type."""
        type_sorted = sorted(dict_, key=self._type_sort_key)
        type_and_value_sorted = []
        for _, group in itertools.groupby(type_sorted, self._type_sort_key):
            type_and_value_sorted.extend(sorted(group))
        return type_and_value_sorted

    def _type_sort_key(self, key):
        return repr(type(key))

    def _format_metadata_key_value(self, key, value):
        """Format metadata key:value, wrapping across lines if necessary."""
        key_fmt = self._format_key(key)

        supported_type = True
        if isinstance(value, (six.text_type, six.binary_type)):
            # for stringy values, there may be u'' or b'' depending on the type
            # of `value` and version of Python. find the starting quote
            # character so that wrapped text will line up with that instead of
            # the string literal prefix character. for example:
            #
            #     'foo': u'abc def ghi
            #              jkl mno'
            value_repr = repr(value)
            extra_indent = 1
            if not (value_repr.startswith("'") or value_repr.startswith('"')):
                extra_indent = 2
        # handles any number, this includes bool
        elif value is None or isinstance(value, numbers.Number):
            value_repr = repr(value)
            extra_indent = 0
        else:
            supported_type = False

        if not supported_type or len(value_repr) > 140:
            value_repr = str(type(value))
            # extra indent of 1 so that wrapped text lines up past the bracket:
            #
            #     'foo': <type
            #             'dict'>
            extra_indent = 1

        return self._wrap_text_with_indent(value_repr, key_fmt, extra_indent)

    def _format_key(self, key):
        """Format metadata key.

        Includes initial indent and trailing colon and space:

            <indent>'foo':<space>

        """
        key_fmt = self._indent + repr(key)
        supported_types = (six.text_type, six.binary_type, numbers.Number,
                           type(None))
        if len(key_fmt) > (self._width / 2) or not isinstance(key,
                                                              supported_types):
            key_fmt = self._indent + str(type(key))
        return '%s: ' % key_fmt

    def _wrap_text_with_indent(self, text, initial_text, extra_indent):
        """Wrap text across lines with an initial indentation.

        For example:

            'foo': 'abc def
                    ghi jkl
                    mno pqr'

        <indent>'foo':<space> is `initial_text`. `extra_indent` is 1. Wrapped
        lines are indented such that they line up with the start of the
        previous line of wrapped text.

        """
        return textwrap.wrap(
            text, width=self._width, expand_tabs=False,
            initial_indent=initial_text,
            subsequent_indent=' ' * (len(initial_text) + extra_indent))

    def _format_positional_metadata_column(self, key, dtype):
        key_fmt = self._format_key(key)
        dtype_fmt = '<dtype: %s>' % str(dtype)
        return self._wrap_text_with_indent(dtype_fmt, key_fmt, 1)

    def _find_optimal_seq_chunking(self):
        """Find the optimal number of sequence chunks to fit on a single line.

        Returns the number of lines the sequence will occupy, the number of
        sequence characters displayed on each line, and the column width
        necessary to display position info using the optimal number of sequence
        chunks.

        """
        # strategy: use an iterative approach to find the optimal number of
        # sequence chunks per line. start with a single chunk and increase
        # until the max line width is exceeded. when this happens, the previous
        # number of chunks is optimal
        num_lines = 0
        num_chars = 0
        column_width = 0

        num_chunks = 1
        not_exceeded = True
        while not_exceeded:
            line_len, new_chunk_info = self._compute_chunked_seq_line_len(
                num_chunks)
            not_exceeded = line_len <= self._width
            if not_exceeded:
                num_lines, num_chars, column_width = new_chunk_info
                num_chunks += 1
        return num_lines, num_chars, column_width

    def _compute_chunked_seq_line_len(self, num_chunks):
        """Compute line length based on a number of chunks."""
        num_chars = num_chunks * self._chunk_size

        # ceil to account for partial line
        num_lines = int(math.ceil(len(self._seq) / num_chars))

        # position column width is fixed width, based on the number of
        # characters necessary to display the position of the final line (all
        # previous positions will be left justified using this width)
        column_width = len('%d ' % ((num_lines - 1) * num_chars))

        # column width + number of sequence characters + spaces between chunks
        line_len = column_width + num_chars + (num_chunks - 1)
        return line_len, (num_lines, num_chars, column_width)

    def _format_chunked_seq(self, line_idxs, num_chars, column_width):
        """Format specified lines of chunked sequence data."""
        lines = []
        for line_idx in line_idxs:
            seq_idx = line_idx * num_chars
            chars = str(self._seq[seq_idx:seq_idx+num_chars])
            chunked_chars = chunk_str(chars, self._chunk_size, ' ')
            lines.append(('%d' % seq_idx).ljust(column_width) + chunked_chars)
        return lines
