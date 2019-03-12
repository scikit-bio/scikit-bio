# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
import collections
import numbers
from contextlib import contextmanager

import numpy as np
import pandas as pd

import skbio.sequence.distance
from skbio._base import SkbioObject
from skbio.metadata._mixin import (MetadataMixin, PositionalMetadataMixin,
                                   IntervalMetadataMixin)
from skbio.metadata import IntervalMetadata
from skbio.sequence._repr import _SequenceReprBuilder
from skbio.util._decorator import (stable, experimental, classonlymethod,
                                   overrides)


class Sequence(MetadataMixin, PositionalMetadataMixin, IntervalMetadataMixin,
               collections.Sequence, SkbioObject):
    """Store generic sequence data and optional associated metadata.

    ``Sequence`` objects do not enforce an alphabet or grammar and are thus the
    most generic objects for storing sequence data. ``Sequence`` objects do not
    necessarily represent biological sequences. For example, ``Sequence`` can
    be used to represent a position in a multiple sequence alignment.
    Subclasses ``DNA``, ``RNA``, and ``Protein`` enforce the IUPAC character
    set [1]_ for, and provide operations specific to, each respective molecule
    type.

    ``Sequence`` objects consist of the underlying sequence data, as well
    as optional metadata and positional metadata. The underlying sequence
    is immutable, while the metdata and positional metadata are mutable.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\\|S1')
        Characters representing the sequence itself.
    metadata : dict, optional
        Arbitrary metadata which applies to the entire sequence. A shallow copy
        of the ``dict`` will be made (see Examples section below for details).
    positional_metadata : pd.DataFrame consumable, optional
        Arbitrary per-character metadata (e.g., sequence read quality
        scores). Must be able to be passed directly to ``pd.DataFrame``
        constructor. Each column of metadata must be the same length as
        `sequence`. A shallow copy of the positional metadata will be made if
        necessary (see Examples section below for details).
    interval_metadata : IntervalMetadata
        Arbitrary metadata which applies to intervals within a sequence to
        store interval features (such as genes, ncRNA on the sequence).
    lowercase : bool or str, optional
        If ``True``, lowercase sequence characters will be converted to
        uppercase characters. If ``False``, no characters will be converted.
        If a str, it will be treated as a key into the positional metadata of
        the object. All lowercase characters will be converted to uppercase,
        and a ``True`` value will be stored in a boolean array in the
        positional metadata under the key.

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
    >>> from skbio.metadata import IntervalMetadata

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

    >>> metadata = {'authors': ['Alice'], 'desc':'seq desc', 'id':'seq-id'}
    >>> positional_metadata = {'exons': [True, True, False, True],
    ...                        'quality': [3, 3, 4, 10]}
    >>> interval_metadata = IntervalMetadata(4)
    >>> interval = interval_metadata.add([(1, 3)], metadata={'gene': 'sagA'})
    >>> seq = Sequence('ACGT', metadata=metadata,
    ...                positional_metadata=positional_metadata,
    ...                interval_metadata=interval_metadata)
    >>> seq
    Sequence
    -----------------------------
    Metadata:
        'authors': <class 'list'>
        'desc': 'seq desc'
        'id': 'seq-id'
    Positional metadata:
        'exons': <dtype: bool>
        'quality': <dtype: int64>
    Interval metadata:
        1 interval feature
    Stats:
        length: 4
    -----------------------------
    0 ACGT

    **Retrieving underlying sequence data:**

    Retrieve underlying sequence:

    >>> seq.values # doctest: +NORMALIZE_WHITESPACE
    array([b'A', b'C', b'G', b'T'],
          dtype='|S1')

    Underlying sequence immutable:

    >>> seq.values = np.array([b'T', b'C', b'G', b'A'], dtype='|S1')
    Traceback (most recent call last):
        ...
    AttributeError: can't set attribute

    >>> seq.values[0] = b'T'
    Traceback (most recent call last):
        ...
    ValueError: assignment destination is read-only

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

    Retrieve interval metadata:

    >>> seq.interval_metadata   # doctest: +ELLIPSIS
    1 interval feature
    ------------------
    Interval(interval_metadata=<...>, bounds=[(1, 3)], \
fuzzy=[(False, False)], metadata={'gene': 'sagA'})

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
    -----------------------------
    Metadata:
        'authors': <class 'list'>
        'desc': 'seq desc'
        'id': 'new-id'
        'pubmed': 12345
    Stats:
        length: 2
    -----------------------------
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
    ...     {'list': [[], [], [], []], 'quality': [3, 3, 4, 10]})
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

    You can also update the interval metadata. Let's re-create a
    ``Sequence`` object with interval metadata at first:

    >>> seq = Sequence('ACGT')
    >>> interval = seq.interval_metadata.add(
    ...     [(1, 3)], metadata={'gene': 'foo'})

    You can update directly on the ``Interval`` object:

    >>> interval  # doctest: +ELLIPSIS
    Interval(interval_metadata=<...>, bounds=[(1, 3)], \
fuzzy=[(False, False)], metadata={'gene': 'foo'})
    >>> interval.bounds = [(0, 2)]
    >>> interval  # doctest: +ELLIPSIS
    Interval(interval_metadata=<...>, bounds=[(0, 2)], \
fuzzy=[(False, False)], metadata={'gene': 'foo'})

    You can also query and obtain the interval features you are
    interested and then modify them:

    >>> intervals = list(seq.interval_metadata.query(metadata={'gene': 'foo'}))
    >>> intervals[0].fuzzy = [(True, False)]
    >>> print(intervals[0])  # doctest: +ELLIPSIS
    Interval(interval_metadata=<...>, bounds=[(0, 2)], \
fuzzy=[(True, False)], metadata={'gene': 'foo'})
    """
    _number_of_extended_ascii_codes = 256
    # ASCII is built such that the difference between uppercase and lowercase
    # is the 6th bit.
    _ascii_invert_case_bit_offset = 32
    _ascii_lowercase_boundary = 90
    default_write_format = 'fasta'
    __hash__ = None

    @property
    @stable(as_of="0.4.0")
    def values(self):
        r"""Array containing underlying sequence characters.

        Notes
        -----
        This property is not writeable.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('AACGA')
        >>> s.values # doctest: +NORMALIZE_WHITESPACE
        array([b'A', b'A', b'C', b'G', b'A'],
              dtype='|S1')

        """
        return self._bytes.view('|S1')

    @property
    def __array_interface__(self):
        r"""Array interface for compatibility with numpy.

        This property allows a ``Sequence`` object to share its underlying data
        buffer (``Sequence.values``) with numpy. See [1]_ for more details.

        References
        ----------
        .. [1] http://docs.scipy.org/doc/numpy/reference/arrays.interface.html

        Examples
        --------
        >>> import numpy as np
        >>> from skbio import Sequence
        >>> seq = Sequence('ABC123')
        >>> np.asarray(seq) # doctest: +NORMALIZE_WHITESPACE
        array([b'A', b'B', b'C', b'1', b'2', b'3'],
              dtype='|S1')

        """
        return self.values.__array_interface__

    @property
    @experimental(as_of="0.4.1")
    def observed_chars(self):
        r"""Set of observed characters in the sequence.

        Notes
        -----
        This property is not writeable.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('AACGAC')
        >>> s.observed_chars == {'G', 'A', 'C'}
        True

        """
        return set(str(self))

    @property
    def _string(self):
        return self._bytes.tostring()

    @classonlymethod
    @experimental(as_of="0.4.1")
    def concat(cls, sequences, how='strict'):
        r"""Concatenate an iterable of ``Sequence`` objects.

        Parameters
        ----------
        sequences : iterable (Sequence)
            An iterable of ``Sequence`` objects or appropriate subclasses.
        how : {'strict', 'inner', 'outer'}, optional
            How to intersect the `positional_metadata` of the sequences.
            If 'strict': the `positional_metadata` must have the exact same
            columns; 'inner': an inner-join of the columns (only the shared set
            of columns are used); 'outer': an outer-join of the columns
            (all columns are used: missing values will be padded with NaN).

        Returns
        -------
        Sequence
            The returned sequence will be an instance of the class which
            called this class-method.

        Raises
        ------
        ValueError
            If `how` is not one of: 'strict', 'inner', or 'outer'.
        ValueError
            If `how` is 'strict' and the `positional_metadata` of each sequence
            does not have the same columns.
        TypeError
            If the sequences cannot be cast as the calling class.

        Notes
        -----
        The sequence-wide metadata (``Sequence.metadata``) is not retained
        during concatenation.

        Sequence objects can be cast to a different type only when the new
        type is an ancestor or child of the original type. Casting between
        sibling types is not allowed, e.g. ``DNA`` -> ``RNA`` is not
        allowed, but ``DNA`` -> ``Sequence`` or ``Sequence`` -> ``DNA``
        would be.

        Examples
        --------
        Concatenate two DNA sequences into a new DNA object:

        >>> from skbio import DNA, Sequence
        >>> s1 = DNA("ACGT")
        >>> s2 = DNA("GGAA")
        >>> DNA.concat([s1, s2])
        DNA
        --------------------------
        Stats:
            length: 8
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 50.00%
        --------------------------
        0 ACGTGGAA

        Concatenate DNA sequences into a Sequence object (type coercion):

        >>> Sequence.concat([s1, s2])
        Sequence
        -------------
        Stats:
            length: 8
        -------------
        0 ACGTGGAA

        Positional metadata is conserved:

        >>> s1 = DNA('AcgT', lowercase='one')
        >>> s2 = DNA('GGaA', lowercase='one',
        ...          positional_metadata={'two': [1, 2, 3, 4]})
        >>> result = DNA.concat([s1, s2], how='outer')
        >>> result
        DNA
        ---------------------------
        Positional metadata:
            'one': <dtype: bool>
            'two': <dtype: float64>
        Stats:
            length: 8
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 50.00%
        ---------------------------
        0 ACGTGGAA
        >>> result.positional_metadata
             one  two
        0  False  NaN
        1   True  NaN
        2   True  NaN
        3  False  NaN
        4  False  1.0
        5  False  2.0
        6   True  3.0
        7  False  4.0

        """
        if how not in {'strict', 'inner', 'outer'}:
            raise ValueError("`how` must be 'strict', 'inner', or 'outer'.")

        seqs = list(sequences)
        if len(seqs) == 0:
            return cls("")

        for seq in seqs:
            seq._assert_can_cast_to(cls)

        if how == 'strict':
            how = 'inner'
            cols = set()
            for s in seqs:
                if s.has_positional_metadata():
                    cols.add(frozenset(s.positional_metadata))
                else:
                    cols.add(frozenset())
            if len(cols) > 1:
                raise ValueError("The positional metadata of the sequences do"
                                 " not have matching columns. Consider setting"
                                 " how='inner' or how='outer'")
        seq_data = []
        pm_data = []
        for seq in seqs:
            seq_data.append(seq._bytes)
            pm_data.append(seq.positional_metadata)
            if not seq.has_positional_metadata():
                del seq.positional_metadata

        pm = pd.concat(pm_data, join=how, ignore_index=True, sort=True)
        bytes_ = np.concatenate(seq_data)

        im = IntervalMetadata.concat(i.interval_metadata for i in seqs)

        return cls(bytes_, positional_metadata=pm, interval_metadata=im)

    @classmethod
    def _assert_can_cast_to(cls, target):
        if not (issubclass(cls, target) or issubclass(target, cls)):
            raise TypeError("Cannot cast %r as %r." %
                            (cls.__name__, target.__name__))

    @overrides(PositionalMetadataMixin)
    def _positional_metadata_axis_len_(self):
        return len(self)

    @overrides(IntervalMetadataMixin)
    def _interval_metadata_axis_len_(self):
        return len(self)

    @stable(as_of="0.4.0")
    def __init__(self, sequence, metadata=None, positional_metadata=None,
                 interval_metadata=None, lowercase=False):
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
            # Sequence casting is acceptable between direct
            # decendants/ancestors
            sequence._assert_can_cast_to(type(self))

            if metadata is None and sequence.has_metadata():
                metadata = sequence.metadata
            if (positional_metadata is None and
                    sequence.has_positional_metadata()):
                positional_metadata = sequence.positional_metadata
            if (interval_metadata is None and
                    sequence.has_interval_metadata()):
                interval_metadata = sequence.interval_metadata
            sequence = sequence._bytes
            self._owns_bytes = False
            self._set_bytes(sequence)
        else:
            # Encode as ascii to raise UnicodeEncodeError if necessary.
            if isinstance(sequence, str):
                sequence = sequence.encode("ascii")
            s = np.frombuffer(sequence, dtype=np.uint8)

            # There are two possibilities (to our knowledge) at this point:
            # Either the sequence we were given was something string-like,
            # (else it would not have made it past frombuffer), or it was a
            # numpy scalar, and so our length must be 1.
            if isinstance(sequence, np.generic) and len(s) != 1:
                raise TypeError("Can cannot create a sequence with %r" %
                                type(sequence).__name__)

            sequence = s
            self._owns_bytes = False

            self._set_bytes(sequence)

        MetadataMixin._init_(self, metadata=metadata)
        PositionalMetadataMixin._init_(
            self, positional_metadata=positional_metadata)
        IntervalMetadataMixin._init_(
            self, interval_metadata=interval_metadata)

        if lowercase is False:
            pass
        elif lowercase is True or isinstance(lowercase, str):
            lowercase_mask = self._bytes > self._ascii_lowercase_boundary
            self._convert_to_uppercase(lowercase_mask)

            # If it isn't True, it must be a string_type
            if not (lowercase is True):
                self.positional_metadata[lowercase] = lowercase_mask
        else:
            raise TypeError("lowercase keyword argument expected a bool or "
                            "string, but got %s" % type(lowercase))

    def _set_bytes_contiguous(self, sequence):
        r"""Munge the sequence data into a numpy array of dtype uint8."""
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

    def _convert_to_uppercase(self, lowercase):
        if np.any(lowercase):
            with self._byte_ownership():
                self._bytes[lowercase] ^= self._ascii_invert_case_bit_offset

    @stable(as_of="0.4.0")
    def __contains__(self, subsequence):
        r"""Determine if a subsequence is contained in this sequence.

        Parameters
        ----------
        subsequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            The putative subsequence.

        Returns
        -------
        bool
            Indicates whether `subsequence` is contained in this sequence.

        Raises
        ------
        TypeError
            If `subsequence` is a ``Sequence`` object with a different type
            than this sequence.

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

    @stable(as_of="0.4.0")
    def __eq__(self, other):
        r"""Determine if this sequence is equal to another.

        Sequences are equal if they are *exactly* the same type and their
        sequence characters, metadata, and positional metadata are the same.

        Parameters
        ----------
        other : Sequence
            Sequence to test for equality against.

        Returns
        -------
        bool
            Indicates whether this sequence is equal to `other`.

        Examples
        --------
        Define two ``Sequence`` objects that have the same underlying sequence
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

        Define another sequence object with a different sequence of characters
        than the previous two sequence objects:

        >>> u = Sequence('ACGA')
        >>> u == t
        False

        Define a sequence with the same sequence of characters as ``u`` but
        with different metadata, positional metadata, and interval metadata:

        >>> v = Sequence('ACGA', metadata={'id': 'abc'},
        ...              positional_metadata={'quality':[1, 5, 3, 3]})
        >>> _ = v.interval_metadata.add([(0, 1)])

        The two sequences are not considered equal because their metadata,
        positional metadata, and interval metadata do not match:

        >>> u == v
        False

        """
        # checks ordered from least to most expensive
        if self.__class__ != other.__class__:
            return False

        if not MetadataMixin._eq_(self, other):
            return False

        if self._string != other._string:
            return False

        if not PositionalMetadataMixin._eq_(self, other):
            return False

        if not IntervalMetadataMixin._eq_(self, other):
            return False

        return True

    @stable(as_of="0.4.0")
    def __ne__(self, other):
        r"""Determine if this sequence is not equal to another.

        Sequences are not equal if they are not *exactly* the same type, or
        their sequence characters, metadata, or positional metadata differ.

        Parameters
        ----------
        other : Sequence
            Sequence to test for inequality against.

        Returns
        -------
        bool
            Indicates whether this sequence is not equal to `other`.

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

    @stable(as_of="0.4.0")
    def __getitem__(self, indexable):
        r"""Slice this sequence.

        Notes
        -----
        This drops the ``self.interval_metadata`` from the returned
        new ``Sequence`` object.

        Parameters
        ----------
        indexable : int, slice, iterable (int and slice), 1D array_like (bool)
            The position(s) to return from this sequence. If `indexable` is an
            iterable of integers, these are assumed to be indices in the
            sequence to keep. If `indexable` is a 1D ``array_like`` of
            booleans, these are assumed to be the positions in the sequence to
            keep.

        Returns
        -------
        Sequence
            New sequence containing the position(s) specified by `indexable` in
            this sequence. Positional metadata will be sliced in the same
            manner and included in the returned sequence. `metadata` is
            included in the returned sequence.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')

        Obtain a single character from the sequence:

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
            ((not isinstance(indexable, str)) and
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
                        positional_metadata = pd.concat(pos_md_slices,
                                                        sort=True)

                    metadata = None
                    if self.has_metadata():
                        metadata = self.metadata

                    return self._constructor(
                        sequence=seq,
                        metadata=metadata,
                        positional_metadata=positional_metadata)
        elif (isinstance(indexable, str) or
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

        metadata = None
        if self.has_metadata():
            metadata = self.metadata

        return self._constructor(
            sequence=seq,
            metadata=metadata,
            positional_metadata=positional_metadata)

    def _slice_positional_metadata(self, indexable):
        if self.has_positional_metadata():
            if _is_single_index(indexable):
                index = _single_index_to_slice(indexable)
            else:
                index = indexable
            return self.positional_metadata.iloc[index]
        else:
            return None

    @stable(as_of="0.4.0")
    def __len__(self):
        r"""Return the number of characters in this sequence.

        Returns
        -------
        int
            The length of this sequence.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> len(s)
        4

        """
        return self._bytes.size

    @stable(as_of="0.4.0")
    def __bool__(self):
        r"""Returns truth value (truthiness) of sequence.

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

    @stable(as_of="0.4.0")
    def __iter__(self):
        r"""Iterate over positions in this sequence.

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

    @stable(as_of="0.4.0")
    def __reversed__(self):
        r"""Iterate over positions in this sequence in reverse order.

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

    @stable(as_of="0.4.0")
    def __str__(self):
        r"""Return sequence characters as a string.

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

    @stable(as_of="0.4.0")
    def __repr__(self):
        r"""Return a string representation of this sequence object.

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
        * interval metadata: the number of interval features will be displayed.
        * sequence stats (e.g., length)
        * up to five lines of chunked sequence data. Each line of chunked
          sequence data displays the current position in the sequence

        Returns
        -------
        str
            String representation of this sequence object.

        Notes
        -----
        Subclasses can override Sequence._repr_stats to provide custom
        statistics.

        Examples
        --------
        Short sequence without metadata:

        >>> from skbio import Sequence
        >>> from skbio.metadata._interval import IntervalMetadata
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

        Sequence with metadata, positional metadata, and interval metadata:

        >>> metadata = {
        ...     'id': 'seq-id',
        ...     'description': 'description of the sequence, wrapping across '
        ...     'lines if it\'s too long',
        ...     'authors': ['Alice', 'Bob', 'Carol'],
        ...     'year': 2015,
        ...     'published': True
        ... }
        >>> positional_metadata = {
        ...     'exons': [True, True, False, True],
        ...     'quality': [3, 10, 11, 10]
        ... }
        >>> seq = Sequence('ACGT', metadata=metadata,
        ...          positional_metadata=positional_metadata)
        >>> _ = seq.interval_metadata.add([(0, 2)], metadata={'gene': 'sagA'})
        >>> seq
        Sequence
        ----------------------------------------------------------------------
        Metadata:
            'authors': <class 'list'>
            'description': "description of the sequence, wrapping across lines
                            if it's too long"
            'id': 'seq-id'
            'published': True
            'year': 2015
        Positional metadata:
            'exons': <dtype: bool>
            'quality': <dtype: int64>
        Interval metadata:
            1 interval feature
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
        r"""Define statistics to display in the sequence's repr.

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

    @stable(as_of="0.4.0")
    def __copy__(self):
        r"""Return a shallow copy of this sequence.

        See Also
        --------
        copy

        Notes
        -----
        This method is equivalent to ``seq.copy(deep=False)``.

        """
        return self._copy(False, {})

    @stable(as_of="0.4.0")
    def __deepcopy__(self, memo):
        r"""Return a deep copy of this sequence.

        See Also
        --------
        copy

        Notes
        -----
        This method is equivalent to ``seq.copy(deep=True)``.

        """
        return self._copy(True, memo)

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
        bytes_ = np.copy(self._bytes)

        seq_copy = self._constructor(sequence=bytes_, metadata=None,
                                     positional_metadata=None,
                                     interval_metadata=None)

        if deep:
            seq_copy._metadata = MetadataMixin._deepcopy_(self, memo)
            seq_copy._positional_metadata = \
                PositionalMetadataMixin._deepcopy_(self, memo)
            seq_copy._interval_metadata = IntervalMetadataMixin._deepcopy_(
                self, memo)
        else:
            seq_copy._metadata = MetadataMixin._copy_(self)
            seq_copy._positional_metadata = \
                PositionalMetadataMixin._copy_(self)
            seq_copy._interval_metadata = IntervalMetadataMixin._copy_(
                self)

        return seq_copy

    @stable(as_of='0.4.0')
    def lowercase(self, lowercase):
        r"""Return a case-sensitive string representation of the sequence.

        Parameters
        ----------
        lowercase: str or boolean vector
            If lowercase is a boolean vector, it is used to set sequence
            characters to lowercase in the output string. True values in the
            boolean vector correspond to lowercase characters. If lowercase
            is a str, it is treated like a key into the positional metadata,
            pointing to a column which must be a boolean vector.
            That boolean vector is then used as described previously.

        Returns
        -------
        str
            String representation of sequence with specified characters set to
            lowercase.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('ACGT')
        >>> s.lowercase([True, True, False, False])
        'acGT'
        >>> s = Sequence('ACGT',
        ...              positional_metadata={
        ...                 'exons': [True, False, False, True]})
        >>> s.lowercase('exons')
        'aCGt'

        Constructor automatically populates a column in positional metadata
        when the ``lowercase`` keyword argument is provided with a column name:

        >>> s = Sequence('ACgt', lowercase='introns')
        >>> s.lowercase('introns')
        'ACgt'
        >>> s = Sequence('ACGT', lowercase='introns')
        >>> s.lowercase('introns')
        'ACGT'

        """
        index = self._munge_to_index_array(lowercase)
        outbytes = self._bytes.copy()
        outbytes[index] ^= self._ascii_invert_case_bit_offset
        return str(outbytes.tostring().decode('ascii'))

    @stable(as_of="0.4.0")
    def count(self, subsequence, start=None, end=None):
        r"""Count occurrences of a subsequence in this sequence.

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
            Number of occurrences of `subsequence` in this sequence.

        Raises
        ------
        ValueError
            If `subsequence` is of length 0.
        TypeError
            If `subsequence` is a ``Sequence`` object with a different type
            than this sequence.

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

    @experimental(as_of="0.5.0")
    def replace(self, where, character):
        r"""Replace values in this sequence with a different character.

        Parameters
        ----------
        where : 1D array_like (bool) or iterable (slices or ints) or str
            Indicates positions in the sequence to replace with `character`.
            Can be a boolean vector, an iterable of indices/slices, or a
            string that is a key in `positional_metadata` pointing to a
            boolean vector.
        character : str or bytes
            Character that will replace chosen items in this sequence.

        Returns
        -------
        Sequence
            Copy of this sequence, with chosen items replaced with chosen
            character. All metadata is retained.

        Examples
        --------
        Let's create and display a Sequence:

        >>> from skbio import Sequence
        >>> sequence = Sequence('GGTACCAACG')
        >>> str(sequence)
        'GGTACCAACG'

        Let's call ``replace`` on the Sequence using a boolean vector for
        ``where`` and assign it to a new variable:

        >>> seq = sequence.replace([False, False, False, True, False, False,
        ...                         True, True, False, False], '-')

        Let's take a look at the new Sequence:

        >>> str(seq)
        'GGT-CC--CG'

        Other types of input are accepted by the ``where`` parameter. Let's
        pass in a list of indices and slices that is equivalent to the boolean
        vector we used previously:

        >>> str(seq) == str(sequence.replace([3, slice(6, 8)], '-'))
        True

        ``where`` also accepts a boolean vector contained in
        ``Sequence.positional_metadata``:

        >>> sequence.positional_metadata = {'where':
        ...                                 [False, False, False, True, False,
        ...                                  False, True, True, False, False]}

        Let's pass in the key ``'where'`` and compare to ``seq``:

        >>> str(seq) == str(sequence.replace('where', '-'))
        True

        """
        if type(character) is not bytes:
            character = character.encode('ascii')
        character = ord(character)
        index = self._munge_to_index_array(where)
        seq_bytes = self._bytes.copy()
        seq_bytes[index] = character

        metadata = None
        if self.has_metadata():
            metadata = self.metadata

        positional_metadata = None
        if self.has_positional_metadata():
            positional_metadata = self.positional_metadata

        interval_metadata = None
        if self.has_interval_metadata():
            interval_metadata = self.interval_metadata
        # Use __class__ instead of _constructor so that validations are
        # performed for subclasses (the user could have introduced invalid
        # characters).
        return self.__class__(seq_bytes, metadata=metadata,
                              positional_metadata=positional_metadata,
                              interval_metadata=interval_metadata)

    @stable(as_of="0.4.0")
    def index(self, subsequence, start=None, end=None):
        r"""Find position where subsequence first occurs in the sequence.

        Parameters
        ----------
        subsequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Subsequence to search for in this sequence.
        start : int, optional
            The position at which to start searching (inclusive).
        end : int, optional
            The position at which to stop searching (exclusive).

        Returns
        -------
        int
            Position where `subsequence` first occurs in this sequence.

        Raises
        ------
        ValueError
            If `subsequence` is not present in this sequence.
        TypeError
            If `subsequence` is a ``Sequence`` object with a different type
            than this sequence.

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

    @experimental(as_of="0.4.0")
    def distance(self, other, metric=None):
        r"""Compute the distance to another sequence.

        Parameters
        ----------
        other : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Sequence to compute the distance to. If `other` is a ``Sequence``
            object, it must be the same type as this sequence. Other input
            types will be converted into a ``Sequence`` object of the same type
            as this sequence.
        metric : function, optional
            Function used to compute the distance between this sequence and
            `other`. If ``None`` (the default), Hamming distance will be used
            (:func:`skbio.sequence.distance.hamming`). `metric` should take two
            ``skbio.Sequence`` objects and return a ``float``. The sequence
            objects passed to `metric` will be the same type as this sequence.
            See :mod:`skbio.sequence.distance` for other predefined metrics
            that can be supplied via `metric`.

        Returns
        -------
        float
            Distance between this sequence and `other` as defined by `metric`.

        Raises
        ------
        TypeError
            If `other` is a ``Sequence`` object with a different type than this
            sequence.

        See Also
        --------
        skbio.sequence.distance
        fraction_diff
        fraction_same

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('AGUC')

        Compute Hamming distance (the default metric):

        >>> s.distance(t)
        0.25

        Use a custom metric:

        >>> def custom_metric(s1, s2): return 0.42
        >>> s.distance(t, custom_metric)
        0.42

        """
        # TODO refactor this method to accept a name (string) of the distance
        # metric to apply and accept **kwargs
        other = self._munge_to_self_type(other, 'distance')
        if metric is None:
            metric = skbio.sequence.distance.hamming
        return float(metric(self, other))

    @stable(as_of="0.4.0")
    def matches(self, other):
        r"""Find positions that match with another sequence.

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
            If `other` is a ``Sequence`` object with a different type than this
            sequence.

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

    @stable(as_of="0.4.0")
    def mismatches(self, other):
        r"""Find positions that do not match with another sequence.

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
            If `other` is a ``Sequence`` object with a different type than this
            sequence.

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

    @stable(as_of="0.4.0")
    def match_frequency(self, other, relative=False):
        r"""Return count of positions that are the same between two sequences.

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
            If `other` is a ``Sequence`` object with a different type than this
            sequence.

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

    @stable(as_of="0.4.0")
    def mismatch_frequency(self, other, relative=False):
        r"""Return count of positions that differ between two sequences.

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
            If `other` is a ``Sequence`` object with a different type than this
            sequence.

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

    @experimental(as_of="0.4.1")
    def frequencies(self, chars=None, relative=False):
        r"""Compute frequencies of characters in the sequence.

        Parameters
        ----------
        chars : str or set of str, optional
            Characters to compute the frequencies of. May be a ``str``
            containing a single character or a ``set`` of single-character
            strings. If ``None``, frequencies will be computed for all
            characters present in the sequence.
        relative : bool, optional
            If ``True``, return the relative frequency of each character
            instead of its count. If `chars` is provided, relative frequencies
            will be computed with respect to the number of characters in the
            sequence, **not** the total count of characters observed in
            `chars`. Thus, the relative frequencies will not necessarily sum to
            1.0 if `chars` is provided.

        Returns
        -------
        dict
            Frequencies of characters in the sequence.

        Raises
        ------
        TypeError
            If `chars` is not a ``str`` or ``set`` of ``str``.
        ValueError
            If `chars` is not a single-character ``str`` or a ``set`` of
            single-character strings.
        ValueError
            If `chars` contains characters outside the allowable range of
            characters in a ``Sequence`` object.

        See Also
        --------
        kmer_frequencies
        iter_kmers

        Notes
        -----
        If the sequence is empty (i.e., length zero), ``relative=True``,
        **and** `chars` is provided, the relative frequency of each specified
        character will be ``np.nan``.

        If `chars` is not provided, this method is equivalent to, but faster
        than, ``seq.kmer_frequencies(k=1)``.

        If `chars` is not provided, it is equivalent to, but faster than,
        passing ``chars=seq.observed_chars``.

        Examples
        --------
        Compute character frequencies of a sequence:

        >>> from pprint import pprint
        >>> from skbio import Sequence
        >>> seq = Sequence('AGAAGACC')
        >>> freqs = seq.frequencies()
        >>> pprint(freqs) # using pprint to display dict in sorted order
        {'A': 4, 'C': 2, 'G': 2}

        Compute relative character frequencies:

        >>> freqs = seq.frequencies(relative=True)
        >>> pprint(freqs)
        {'A': 0.5, 'C': 0.25, 'G': 0.25}

        Compute relative frequencies of characters A, C, and T:

        >>> freqs = seq.frequencies(chars={'A', 'C', 'T'}, relative=True)
        >>> pprint(freqs)
        {'A': 0.5, 'C': 0.25, 'T': 0.0}

        Note that since character T is not in the sequence we receive a
        relative frequency of 0.0. The relative frequencies of A and C are
        relative to the number of characters in the sequence (8), **not** the
        number of A and C characters (4 + 2 = 6).

        """
        freqs = np.bincount(self._bytes,
                            minlength=self._number_of_extended_ascii_codes)

        if chars is not None:
            chars, indices = self._chars_to_indices(chars)
        else:
            indices, = np.nonzero(freqs)
            # Downcast from int64 to uint8 then convert to str. This is safe
            # because we are guaranteed to have indices in the range 0 to 255
            # inclusive.
            chars = indices.astype(np.uint8).tostring().decode('ascii')

        obs_counts = freqs[indices]
        if relative:
            obs_counts = obs_counts / len(self)

        # Use tolist() for minor performance gain.
        return dict(zip(chars, obs_counts.tolist()))

    def _chars_to_indices(self, chars):
        r"""Helper for Sequence.frequencies."""
        if isinstance(chars, (str, bytes)):
            chars = set([chars])
        elif not isinstance(chars, set):
            raise TypeError(
                "`chars` must be of type `set`, not %r" % type(chars).__name__)

        # Impose an (arbitrary) ordering to `chars` so that we can return
        # `indices` in that same order.
        chars = list(chars)
        indices = []
        for char in chars:
            if not isinstance(char, (str, bytes)):
                raise TypeError(
                    "Each element of `chars` must be string-like, not %r" %
                    type(char).__name__)
            if len(char) != 1:
                raise ValueError(
                    "Each element of `chars` must contain a single "
                    "character (found %d characters)" % len(char))

            index = ord(char)
            if index >= self._number_of_extended_ascii_codes:
                raise ValueError(
                    "Character %r in `chars` is outside the range of "
                    "allowable characters in a `Sequence` object." % char)
            indices.append(index)
        return chars, indices

    @stable(as_of="0.4.0")
    def iter_kmers(self, k, overlap=True):
        r"""Generate kmers of length `k` from this sequence.

        Parameters
        ----------
        k : int
            The kmer length.
        overlap : bool, optional
            Defines whether the kmers should be overlapping or not.

        Yields
        ------
        Sequence
            kmer of length `k` contained in this sequence.

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

        if len(self) == 0 or self.has_positional_metadata():
            # Slower path when sequence is empty or positional metadata needs
            # to be sliced.
            for i in range(0, len(self) - k + 1, step):
                yield self[i:i+k]
        else:
            # Optimized path when positional metadata doesn't need slicing.
            kmers = np.lib.stride_tricks.as_strided(
                self._bytes, shape=(k, count), strides=(1, step)).T

            metadata = None
            if self.has_metadata():
                metadata = self.metadata

            for s in kmers:
                yield self._constructor(
                    sequence=s,
                    metadata=metadata,
                    positional_metadata=None)

    @stable(as_of="0.4.0")
    def kmer_frequencies(self, k, overlap=True, relative=False):
        r"""Return counts of words of length `k` from this sequence.

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
        dict
            Frequencies of words of length `k` contained in this sequence.

        Raises
        ------
        ValueError
            If `k` is less than 1.

        Examples
        --------
        >>> from pprint import pprint
        >>> from skbio import Sequence
        >>> s = Sequence('ACACATTTATTA')
        >>> freqs = s.kmer_frequencies(3, overlap=False)
        >>> pprint(freqs) # using pprint to display dict in sorted order
        {'ACA': 1, 'CAT': 1, 'TTA': 2}
        >>> freqs = s.kmer_frequencies(3, relative=True, overlap=False)
        >>> pprint(freqs)
        {'ACA': 0.25, 'CAT': 0.25, 'TTA': 0.5}

        """
        kmers = self.iter_kmers(k, overlap=overlap)
        freqs = dict(collections.Counter((str(seq) for seq in kmers)))

        if relative:
            if overlap:
                num_kmers = len(self) - k + 1
            else:
                num_kmers = len(self) // k

            relative_freqs = {}
            for kmer, count in freqs.items():
                relative_freqs[kmer] = count / num_kmers
            freqs = relative_freqs

        return freqs

    @stable(as_of="0.4.0")
    def find_with_regex(self, regex, ignore=None):
        r"""Generate slices for patterns matched by a regular expression.

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
        if isinstance(regex, str):
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

    @stable(as_of="0.4.0")
    def iter_contiguous(self, included, min_length=1, invert=False):
        r"""Yield contiguous subsequences based on `included`.

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

    def _constructor(self, **kwargs):
        return self.__class__(**kwargs)

    def _munge_to_index_array(self, sliceable):
        r"""Return an index array from something isomorphic to a boolean vector.

        """
        if isinstance(sliceable, str):
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

    def _munge_to_self_type(self, other, method):
        if isinstance(other, Sequence):
            if type(other) != type(self):
                raise TypeError("Cannot use %s and %s together with `%s`" %
                                (self.__class__.__name__,
                                 other.__class__.__name__, method))
            else:
                return other

        return self.__class__(other)

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
        elif isinstance(other, str):
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
