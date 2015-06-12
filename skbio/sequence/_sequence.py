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
from six import string_types, text_type

import re
import collections
import numbers
from contextlib import contextmanager

import numpy as np
from scipy.spatial.distance import hamming

import pandas as pd

from skbio._base import SkbioObject
from skbio.util._misc import reprnator


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


def _dataframe_with_reset_index(dataframe):
    # By default Pandas adds a column containing the original indices.
    # For now we're just throwing it away using the following.
    return dataframe.reset_index(drop=True)


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
        Arbitrary metadata which applies to the entire sequence.
    positional_metadata : Pandas DataFrame consumable, optional
        Arbitrary per-character metadata. For example, quality data from
        sequencing reads. Must be able to be passed directly to the Pandas
        DataFrame constructor.

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

    Notes
    -----
    ``Sequence`` objects are immutable.

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    Examples
    --------
    >>> from skbio import Sequence
    >>> s = Sequence('GGUCGUGAAGGA')
    >>> s # doctest: +NORMALIZE_WHITESPACE
    Sequence('GGUCGUGAAGGA', length=12, has_metadata=False,
             has_positional_metadata=False)
    >>> t = Sequence('CAT', metadata={'id':'seq-id', 'desc':'seq desc'},
    ...              positional_metadata={'qual':[42, 42, 1]})
    >>> t # doctest: +NORMALIZE_WHITESPACE
    Sequence('CAT', length=3, has_metadata=True,
             has_positional_metadata=True)

    """
    default_write_format = 'fasta'
    __hash__ = None  # TODO revisit hashability when all properties are present

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
        """dict containing metadata which applies to the entire sequence.

        Notes
        -----
        This property is writeable.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('ACGTACGTACGTACGT', metadata={'id': 'seq-id',
        ...                                 'description': 'seq description'})
        >>> s # doctest: +NORMALIZE_WHITESPACE
        Sequence('ACGTACGTACGTACGT', length=16, has_metadata=True,
                 has_positional_metadata=False)
        >>> s.metadata
        {'id': 'seq-id', 'description': 'seq description'}

        """
        return self._metadata

    @property
    def positional_metadata(self):
        """Pandas DataFrame containing metadata on a per-character basis.

        Notes
        -----
        This property is writeable.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('ACGTACGTACGTACGT',
        ...              positional_metadata={'quality': range(16)})
        >>> s # doctest: +NORMALIZE_WHITESPACE
        Sequence('ACGTACGTACGTACGT', length=16, has_metadata=False,
                 has_positional_metadata=True)
        >>> s.positional_metadata['quality'].values
        array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15])

        """
        return self._positional_metadata

    @property
    def _string(self):
        return self._bytes.tostring()

    def __init__(self, sequence, metadata=None,
                 positional_metadata=None):
        if isinstance(sequence, Sequence):
            if metadata is None:
                metadata = sequence.metadata
            if positional_metadata is None:
                positional_metadata = sequence.positional_metadata
            sequence = sequence._bytes

        self._set_sequence(sequence)
        self._set_metadata(metadata)
        self._set_positional_metadata(positional_metadata)

    def _set_sequence(self, sequence):
        """Munge the sequence data into a numpy array of dtype uint8."""
        is_ndarray = isinstance(sequence, np.ndarray)
        if is_ndarray:
            if np.issubdtype(sequence.dtype, np.uint8):
                pass
            elif sequence.dtype == '|S1':
                sequence = sequence.view(np.uint8)
            else:
                raise TypeError(
                    "Can only create sequence from numpy.ndarray of dtype "
                    "np.uint8 or '|S1'. Invalid dtype: %s" %
                    sequence.dtype)

            # numpy doesn't support views of non-contiguous arrays. Since we're
            # making heavy use of views internally, and users may also supply
            # us with a view, make sure we *always* store a contiguous array to
            # avoid hard-to-track bugs. See
            # https://github.com/numpy/numpy/issues/5716
            potential_copy = np.ascontiguousarray(sequence)
            if potential_copy is not sequence:
                self._owns_bytes = True
            else:
                self._owns_bytes = False
            sequence = potential_copy
        else:
            # Python 3 will not raise a UnicodeEncodeError so we force it by
            # encoding it as ascii
            if isinstance(sequence, text_type):
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

        sequence.flags.writeable = False
        self._bytes = sequence

    def _set_metadata(self, metadata):
        if metadata is None:
            metadata = {}
        elif not isinstance(metadata, dict):
            raise TypeError("metadata must be a {}".format(type(dict())))
        self._metadata = metadata

    def _set_positional_metadata(self, positional_metadata):
        if positional_metadata is None:
            # ensure dataframe of proper length
            positional_metadata = pd.DataFrame(index=range(len(self)))

        try:
            self._positional_metadata = pd.DataFrame(positional_metadata)
        except pd.core.common.PandasError as e:
            raise TypeError('Positional metadata invalid. Must be consumable '
                            'by pandas.DataFrame. Original Pandas error '
                            'message: "%s"' % e.message)

        num_rows = len(self.positional_metadata.index)
        if num_rows != len(self):
            raise ValueError(
                "Number of positional metadata values (%d) must match the "
                "number of characters in the sequence (%d)." %
                (num_rows, len(self)))

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

        See Also
        --------
        equals

        Notes
        -----
        See ``Sequence.equals`` for more control over equality testing.

        This method is equivalent to ``self.equals(other)``.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> t = Sequence('GGUCGUGAAGGA')
        >>> s == t
        True
        >>> u = Sequence('GGUCGUGACCGA')
        >>> u == t
        False

        Note that because the positional metadata does not match between ``u``
        and ``v``, they are not considered equal:

        >>> v = Sequence('GGUCGUGACCGA',
        ...              positional_metadata={'quality':[1, 5, 3, 3, 2, 42,
        ...                                              100, 9, 10, 0, 42,
        ...                                              42]})
        >>> u == v
        False

        """
        return self.equals(other)

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

        See Also
        --------
        equals

        Notes
        -----
        See ``Sequence.equals`` for more control over equality testing.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> t = Sequence('GGUCGUGAAGGA')
        >>> s != t
        False
        >>> u = Sequence('GGUCGUGACCGA')
        >>> u != t
        True
        >>> v = Sequence('GGUCGUGACCGA', metadata={'id':'v'})
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

        >>> s[1] # doctest: +NORMALIZE_WHITESPACE
        Sequence('G', length=1, has_metadata=False,
                 has_positional_metadata=False)

        Obtain a slice:

        >>> s[7:] # doctest: +NORMALIZE_WHITESPACE
        Sequence('AAGGA', length=5, has_metadata=False,
                 has_positional_metadata=False)

        Obtain characters at the following indices:

        >>> s[[3, 4, 7, 0, 3]] # doctest: +NORMALIZE_WHITESPACE
        Sequence('CGAGC', length=5, has_metadata=False,
                 has_positional_metadata=False)

        Obtain characters at positions evaluating to `True`:

        >>> s = Sequence('GGUCG')
        >>> index = [True, False, True, 'a' is 'a', False]
        >>> s[index] # doctest: +NORMALIZE_WHITESPACE
        Sequence('GUC', length=3, has_metadata=False,
                 has_positional_metadata=False)

        """
        metadata = self.metadata
        if (not isinstance(indexable, np.ndarray) and
            ((not isinstance(indexable, string_types)) and
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
                    pos_md_slices = list(_slices_from_iter(
                                         self.positional_metadata, index))
                    positional_metadata = \
                        _dataframe_with_reset_index(pd.concat(pos_md_slices))

                    return self._to(sequence=seq, metadata=metadata,
                                    positional_metadata=positional_metadata)
        elif (isinstance(indexable, string_types) or
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
        positional_metadata = None
        if self.has_positional_metadata():
            positional_metadata = self._slice_positional_metadata(indexable)

        return self._to(sequence=seq, metadata=metadata,
                        positional_metadata=positional_metadata)

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
        return bool(self.metadata)

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
        return len(self.positional_metadata.columns) > 0

    def _slice_positional_metadata(self, indexable):
        if _is_single_index(indexable):
            index = _single_index_to_slice(indexable)
        else:
            index = indexable
        return _dataframe_with_reset_index(
            self.positional_metadata.iloc[index])

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
        """Return a string representation of the biological sequence object.

        Returns
        -------
        str
            String representation of the biological sequence object.

        Notes
        -----
        String representation contains the class name, the first six characters
        of the sequence, followed by ellipses, followed by the last six
        characters of the sequence (or the full sequence without
        ellipses, if the sequence is less than 21 characters long), followed by
        the sequence length, followed by flags indicating whether the sequence
        has metadata and/or positional metadata.

        Examples
        --------
        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> s # doctest: +NORMALIZE_WHITESPACE
        Sequence('GGUCGUGAAGGA', length=12, has_metadata=False,
                 has_positional_metadata=False)
        >>> t = Sequence('ACGT')
        >>> t # doctest: +NORMALIZE_WHITESPACE
        Sequence('ACGT', length=4, has_metadata=False,
                 has_positional_metadata=False)
        >>> t # doctest: +NORMALIZE_WHITESPACE
        Sequence('ACGT', length=4, has_metadata=False,
                 has_positional_metadata=False)
        >>> Sequence('GGUCGUGAAAAAAAAAAAAGGA') # doctest: +NORMALIZE_WHITESPACE
        Sequence('GGUCGU ... AAAGGA', length=22, has_metadata=False,
                 has_positional_metadata=False)
        >>> Sequence('ACGT',
        ...          metadata={id:'seq1'}) # doctest: +NORMALIZE_WHITESPACE
        Sequence('ACGT', length=4, has_metadata=True,
                 has_positional_metadata=False)
        """

        start = self.__class__.__name__ + "("
        end = ")"

        tokens = []

        tokens.append(self._format_str(self))
        tokens.append("length=%d" % len(self))
        tokens.append("has_metadata=%s" % self.has_metadata())
        tokens.append("has_positional_metadata=%s" %
                      self.has_positional_metadata())

        return reprnator(start, tokens, end)

    def _format_str(self, s):
        s = repr(str(s))
        if len(s) > 20:
            return "%s ... %s" % (s[:7], s[-7:])
        return s

    def _format_list(self, l):
        l = list(l)
        if len(l) > 13:
            return "[%s, ..., %s]" % (repr(l[:6])[1:-1], repr(l[-6:])[1:-1])
        return "%r" % l

    def equals(self, other, ignore=None):
        """Determine if the biological sequence is equal to another.

        By default, biological sequences are equal if they are *exactly* the
        same type and their sequence characters, metadata, and positional
        metadata are the same.

        Parameters
        ----------
        other : Sequence
            Sequence to test for equality against.
        ignore : iterable (str), optional
            List of features to ignore in the equality test. By default, all
            features must be the same for two biological sequences to be
            considered equal. Features that can be ignored are ``'type'``,
            ``'sequence'``, ``'metadata'``, and ``'positional_metadata'``.

        Returns
        -------
        bool
            Indicates whether the biological sequence is equal to `other`.

        Examples
        --------
        Define two biological sequences that have the same underlying sequence
        of characters:

        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> t = Sequence('GGUCGUGAAGGA')

        The two sequences are considered equal because they are the same type,
        their underlying sequence of characters are the same, and their
        optional attributes (metadata and positional_metadata) were not
        provided:

        >>> s.equals(t)
        True
        >>> t.equals(s)
        True

        Define another biological sequence with a different sequence of
        characters than the previous two biological sequences:

        >>> u = Sequence('GGUCGUGACCGA')
        >>> u.equals(t)
        False

        Define a biological sequence with the same sequence of characters as
        ``u``, but with different metadata and positional metadata:

        >>> v = Sequence('GGUCGUGACCGA', metadata={'id':'abc'},
        ...               positional_metadata={'quality':[1, 5, 3, 3, 2, 42,
        ...                                               100, 9, 10, 55, 42,
        ...                                               42]})

        By default, the two sequences are *not* considered equal because their
        metadata and positional metadata do not match:

        >>> u.equals(v)
        False

        By specifying that the metadata and positonal metadata should be
        ignored, they now compare equal:

        >>> u.equals(v, ignore=['metadata', 'positional_metadata'])
        True

        """
        if ignore is None:
            ignore = {}

        # Checks are ordered from least to most expensive.
        if 'type' not in ignore and self.__class__ != other.__class__:
            return False

        if 'sequence' not in ignore and self._string != other._string:
            return False

        if 'metadata' not in ignore and self.metadata != other.metadata:
            return False

        if ('positional_metadata' not in ignore and not
                self.positional_metadata.equals(other.positional_metadata)):
            return False

        return True

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

        step = 1 if overlap else k

        for i in range(0, len(self) - k + 1, step):
            yield self[i:i+k]

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
        if isinstance(regex, string_types):
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

    def _to(self, **kwargs):
        """Return a copy of the current biological sequence.

        Returns a copy of the current biological sequence, optionally with
        updated attributes specified as keyword arguments.

        Parameters
        ----------
        kwargs : dict, optional
            Keyword arguments passed to the ``Sequence`` (or
            subclass) constructor. The returned copy will have its attributes
            updated based on the values in `kwargs`. If an attribute is
            missing, the copy will keep the same attribute as the current
            biological sequence. Valid attribute names are `'sequence'`,
            `'metadata'`, `'positional_metadata'`. Default behavior is to
            return a copy of the current biological sequence without changing
            any attributes.

        Returns
        -------
        Sequence
            Copy of the current biological sequence, optionally with updated
            attributes based on `kwargs`. Will be the same type as the current
            biological sequence (`self`).

        Notes
        -----
        This is a shallow copy, but since biological sequences are immutable,
        it is conceptually the same as a deep copy.

        This method is the preferred way of creating new instances from an
        existing biological sequence, instead of calling
        ``self.__class__(...)``, as the latter can be error-prone (e.g.,
        it's easy to forget to propagate attributes to the new instance).

        Examples
        --------
        Create a biological sequence:

        >>> from skbio import Sequence
        >>> seq = Sequence('AACCGGTT',
        ...                metadata={'id':'id1'},
        ...                positional_metadata={
        ...                    'qual':[4, 2, 22, 23, 1, 1, 1, 9]
        ...                })

        Create a copy of ``seq``, keeping the same underlying sequence of
        characters and quality scores, while updating the metadata:

        >>> new_seq = seq._to(metadata={'id':'new-id'})

        Note that the copied biological sequence's underlying sequence and
        positional metadata are the same as ``seq``:

        >>> str(new_seq)
        'AACCGGTT'
        >>> new_seq.positional_metadata['qual'].values
        array([ 4,  2, 22, 23,  1,  1,  1,  9])

        The metadata has been updated:

        >>> new_seq.metadata['id']
        'new-id'

        The original biological sequence's metadata has not been changed:

        >>> seq.metadata['id']
        'id1'

        """
        defaults = {'sequence': self._bytes,
                    'metadata': self.metadata,
                    'positional_metadata': self._positional_metadata}
        defaults.update(kwargs)
        return self._constructor(**defaults)

    def _constructor(self, **kwargs):
        return self.__class__(**kwargs)

    def _munge_to_index_array(self, sliceable):
        """Return an index array from something isomorphic to a boolean vector.

        """
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
        if isinstance(other, string_types):
            return six.b(other)
        return self._munge_to_sequence(other, method)._string

    @contextmanager
    def _byte_ownership(self):
        if not self._owns_bytes:
            self._bytes = self._bytes.copy()
            self._owns_bytes = True

        self._bytes.flags.writeable = True
        yield
        self._bytes.flags.writeable = False
