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
from future.standard_library import hooks
from six import string_types

import collections
import numbers
from contextlib import contextmanager

import numpy as np
from scipy.spatial.distance import hamming

from skbio._base import SkbioObject
from skbio.util._misc import reprnator

with hooks():
    from itertools import zip_longest


class Sequence(collections.Sequence, SkbioObject):
    """Base class for biological sequences.

    Parameters
    ----------
    sequence : python collections.Sequence (e.g., str, list or tuple)
        The biological sequence.
    id : str, optional
        The sequence id (e.g., an accession number).
    description : str, optional
        A description or comment about the sequence (e.g., "green
        fluorescent protein").
    quality : 1-D array_like, int, optional
        Phred quality scores stored as nonnegative integers, one per sequence
        character. If provided, must be the same length as the biological
        sequence. Can be a 1-D ``numpy.ndarray`` of integers, or a structure
        that can be converted to this representation using ``numpy.asarray``. A
        copy will *not* be made if `quality` is already a 1-D ``numpy.ndarray``
        with an ``int`` ``dtype``. The array will be made read-only (i.e., its
        ``WRITEABLE`` flag will be set to ``False``).

    Attributes
    ----------
    sequence
    id
    description
    quality

    Raises
    ------
    skbio.sequence.SequenceError
        If `quality` is not the correct shape.

    See Also
    --------
    NucleotideSequence
    DNA
    RNA

    Notes
    -----
    `Sequence` objects are immutable. Where applicable, methods
    return a new object of the same class.
    Subclasses are typically defined by methods relevant to only a specific
    type of biological sequence, and by containing characters only contained in
    the IUPAC standard character set [1]_ for that molecule type.

    Examples
    --------
    >>> from skbio.sequence import Sequence
    >>> s = Sequence('GGUCGUGAAGGA')
    >>> t = Sequence('GGUCCUGAAGGU')

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    """
    default_write_format = 'fasta'
    __hash__ = None  # TODO revisit hashability when all properties are present

    def __init__(self, sequence, id="", description="", quality=None):
        """4 types to rule them all: char vector, byte vector, Sequence, or
        string_types"""
        if isinstance(sequence, Sequence):
            if id == "":
                id = sequence.id
            if description == "":
                description = sequence.description
            if quality is None:
                quality = sequence.quality

            sequence = sequence._bytes

        self._set_id(id)
        self._set_description(description)
        self._set_sequence(sequence)
        self._set_quality(quality)

    @contextmanager
    def _byte_ownership(self):
        if not self._owns_bytes:
            self._bytes = self._bytes.copy()
            self._owns_bytes = True

        self._bytes.flags.writeable = True
        yield
        self._bytes.flags.writeable = False
        self._refresh_chars()

    def _refresh_chars(self):
        if self._bytes.size == 0:
            self.__chars = self._bytes.view('|S1')
        else:
            self.__chars = self._bytes.view('|S%d' % self._bytes.size)

    @property
    def _string(self):
        if self.__chars.size == 0:
            return ''
        else:
            return self.__chars[0]

    def _set_id(self, id_):
        if isinstance(id_, string_types):
            self._id = id_
        else:
            raise TypeError("Sequence ID %r must be a string type." % (id_,))

    def _set_description(self, description):
        if isinstance(description, string_types):
            self._description = description
        else:
            raise TypeError("Sequence description %r must be a string type." %
                            (description,))

    def _set_sequence(self, sequence):
        """Munge the sequence data into a numpy array."""
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
            sequence = np.fromstring(sequence, dtype=np.uint8)
            self._owns_bytes = True

        sequence.flags.writeable = False

        self._bytes = sequence
        self._refresh_chars()

    def _set_quality(self, quality):
        if quality is not None:
            quality = np.asarray(quality)

            if quality.ndim == 0:
                # We have something scalar-like, so create a single-element
                # vector to store it.
                quality = np.reshape(quality, 1)

            if quality.shape == (0,):
                # cannot safe cast an empty vector from float to int
                cast_type = 'unsafe'
            else:
                cast_type = 'safe'

            quality = quality.astype(int, casting=cast_type, copy=False)
            quality.flags.writeable = False

            if quality.ndim != 1:
                raise ValueError(
                    "Quality scores have %d dimension(s). Quality scores must "
                    "be 1-D." % quality.ndim)

            if len(quality) != len(self):
                raise ValueError(
                    "Number of quality scores (%d) must match the "
                    "number of characters in the sequence (%d)." %
                    (len(quality), len(self)))

            if (quality < 0).any():
                raise ValueError(
                    "Quality scores must be greater than or equal to "
                    "zero.")

        self._quality = quality

    def __contains__(self, other):
        """The in operator.

        Parameters
        ----------
        other : str
            The putative subsequence.

        Returns
        -------
        bool
            Indicates whether `other` is contained in `self`.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> 'GGU' in s
        True
        >>> 'CCC' in s
        False

        .. shownumpydoc

        """
        if isinstance(other, string_types):
            return other in self._string

        if isinstance(other, Sequence) and type(other) != type(self):
            klass = self.__class__.__name__
            oklass = other.__class__.__name__
            raise TypeError("'in <%s>' requires string, numpy array, or %s as"
                            " left operand, not %s" % (klass, klass, oklass))

        return Sequence(other)._string in self._string

    def __eq__(self, other):
        """The equality operator.

        Biological sequences are equal if their sequence is the same and they
        are the same type. Identifier, description, and quality scores
        **are ignored**.

        Parameters
        ----------
        other : `Sequence`
            The sequence to test for equality against.

        Returns
        -------
        bool
            Indicates whether `self` and `other` are equal.

        See Also
        --------
        __ne__
        equals

        Notes
        -----
        See ``Sequence.equals`` for more fine-grained control of
        equality testing.

        This method is equivalent to
        ``self.equals(other, ignore=['id', 'description', 'quality'])``.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> t = Sequence('GGUCGUGAAGGA')
        >>> s == t
        True
        >>> u = Sequence('GGUCGUGACCGA')
        >>> u == t
        False

        Note that even though the quality scores do not match between ``u`` and
        ``v``, they are considered equal:

        >>> v = Sequence('GGUCGUGACCGA',
        ...                        quality=[1, 5, 3, 3, 2, 42, 100, 9, 10, 55,
        ...                                 42, 42])
        >>> u == v
        True

        .. shownumpydoc

        """
        return self.equals(other)

    def __getitem__(self, indexable):
        """The indexing operator.

        Parameters
        ----------
        indexable : int, slice, or sequence of ints
            The position(s) to return from the `Sequence`. If `i` is
            a sequence of ints, these are assumed to be indices in the sequence
            to keep.

        Returns
        -------
        Sequence
            New biological sequence containing the character(s) at position(s)
            `i` in the current `Sequence`. If quality scores are
            present, the quality score at position(s) `i` will be included in
            the returned sequence. ID and description are also included.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')

        Obtain a single character from the biological sequence:

        >>> s[1]
        <Sequence: G (length: 1)>

        Obtain a slice:

        >>> s[7:]
        <Sequence: AAGGA (length: 5)>

        Obtain characters at the following indices:

        >>> s[[3, 4, 7, 0, 3]]
        <Sequence: CGAGC (length: 5)>

        .. shownumpydoc

        """
        qual = None
        if (not isinstance(indexable, np.ndarray) and
            ((not isinstance(indexable, string_types)) and
             hasattr(indexable, '__iter__'))):
            indexable_ = indexable
            indexable = np.asarray(indexable)
            if indexable.dtype == object:
                indexable = list(indexable_)  # TODO: Don't blow out memory
                seq = np.concatenate(list(self._slices_from_iter(self._bytes,
                                                                 indexable)))
                if self._has_quality():
                    qual = np.concatenate(list(self._slices_from_iter(
                        self.quality, indexable)))

                return self.to(sequence=seq, quality=qual)
        elif isinstance(indexable, string_types) or \
                isinstance(indexable, bool):
            raise IndexError("Cannot index with that type: %r" % indexable)

        if (isinstance(indexable, np.ndarray) and
            indexable.dtype == bool and
            len(indexable) != len(self)):
            raise IndexError("An boolean vector index must be the same length"
                            " as the sequence (%d, not %d)." % (len(self),
                                                                len(indexable)))

        seq = self._bytes[indexable]
        if self._has_quality():
            qual = self.quality[indexable]

        return self.to(sequence=seq, quality=qual)

    def _slices_from_iter(self, array, indexables):
        for i in indexables:
            if isinstance(i, slice):
                pass
            elif isinstance(i, numbers.Integral) and not isinstance(i, bool):
                if i == -1:
                    i = slice(i, None)
                else:
                    i = slice(i, i+1)
            else:
                raise IndexError("Cannot slice sequence from iterable "
                                 "containing %r." % i)

            yield array[i]

    def __iter__(self):
        """The iter operator.

        Returns
        -------
        iterator
            Position iterator for the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUC')
        >>> for c in s: print(c)
        G
        G
        U
        C

        .. shownumpydoc

        """
        if self._has_quality():
            qual = self.quality
        else:
            qual = []

        for c, q in zip_longest(self._string, qual, fillvalue=None):
            yield self.to(sequence=c, quality=q)

    def __len__(self):
        """The len operator.

        Returns
        -------
        int
            The length of the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUC')
        >>> len(s)
        4

        .. shownumpydoc

        """
        return self._bytes.size

    def __ne__(self, other):
        """The inequality operator.

        Biological sequences are not equal if their sequence is different or
        they are not the same type. Identifier, description, and quality scores
        **are ignored**.

        Parameters
        ----------
        other : `Sequence`
            The sequence to test for inequality against.

        Returns
        -------
        bool
            Indicates whether `self` and `other` are not equal.

        See Also
        --------
        __eq__
        equals

        Notes
        -----
        See ``Sequence.equals`` for more fine-grained control of
        equality testing.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> t = Sequence('GGUCGUGAAGGA')
        >>> s != t
        False
        >>> u = Sequence('GGUCGUGACCGA')
        >>> u != t
        True

        .. shownumpydoc

        """
        return not (self == other)

    def __repr__(self):
        """The repr method.

        Returns
        -------
        str
            Returns a string representation of the object.

        Notes
        -----
        String representation contains the class name, the first ten characters
        of the sequence followed by ellipses (or the full sequence
        and no ellipses, if the sequence is less than 11 characters long),
        followed by the sequence length.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> repr(s)
        '<Sequence: GGUCGUGAAG... (length: 12)>'
        >>> t = Sequence('ACGT')
        >>> repr(t)
        '<Sequence: ACGT (length: 4)>'
        >>> t
        <Sequence: ACGT (length: 4)>

        .. shownumpydoc

        """
        start = self.__class__.__name__ + "("
        end = ")[0:%d]" % len(self)

        tokens = []

        tokens.append(self._format_str(self))
        if self.id:
            tokens.append("id=" + self._format_str(self.id))
        if self.description:
            tokens.append("description=" + self._format_str(self.description))
        if self._has_quality():
            tokens.append("quality=" + self._format_list(self.quality))

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

    def __reversed__(self):
        """The reversed operator.

        Returns
        -------
        iterator
            Reverse position iterator for the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUC')
        >>> for c in reversed(s): print(c)
        C
        U
        G
        G

        .. shownumpydoc

        """
        return iter(self[::-1])

    def __str__(self):
        """Document me?"""
        return str(self._string)

    @property
    def sequence(self):
        """String containing underlying biological sequence characters.

        A string representing the characters of the biological sequence.

        Notes
        -----
        This property is not writeable.

        """
        return self._bytes.view('|S1')

    @property
    def id(self):
        """ID of the biological sequence.

        A string representing the identifier (ID) of the biological sequence.

        Notes
        -----
        This property is not writeable.

        """
        return self._id

    @property
    def description(self):
        """Description of the biological sequence.

        A string representing the description of the biological sequence.

        Notes
        -----
        This property is not writeable.

        """
        return self._description

    @property
    def quality(self):
        """Quality scores of the characters in the biological sequence.

        A 1-D ``numpy.ndarray`` of nonnegative integers representing Phred
        quality scores for each character in the biological sequence, or
        ``None`` if quality scores are not present.

        Notes
        -----
        This property is not writeable. A copy of the array is *not* returned.
        The array is read-only (i.e., its ``WRITEABLE`` flag is set to
        ``False``).

        """
        return self._quality

    def _has_quality(self):
        """Return bool indicating presence of quality scores in the sequence.

        Returns
        -------
        bool
            ``True`` if the biological sequence has quality scores, ``False``
            otherwise.

        See Also
        --------
        quality

        """
        return self.quality is not None

    def to(self, **kwargs):
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
            `'id'`, `'description'`, and `'quality'`. Default behavior is to
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
        forgetting to propagate attributes to the new instance).

        Examples
        --------
        Create a biological sequence:

        >>> from skbio import Sequence
        >>> seq = Sequence('AACCGGTT', id='id1',
        ...                          description='biological sequence',
        ...                          quality=[4, 2, 22, 23, 1, 1, 1, 9])

        Create a copy of ``seq``, keeping the same underlying sequence of
        characters and quality scores, while updating ID and description:

        >>> new_seq = seq.copy(id='new-id', description='new description')

        Note that the copied biological sequence's underlying sequence and
        quality scores are the same as ``seq``:

        >>> new_seq.sequence
        'AACCGGTT'
        >>> new_seq.quality
        array([ 4,  2, 22, 23,  1,  1,  1,  9])

        The ID and description have been updated:

        >>> new_seq.id
        'new-id'
        >>> new_seq.description
        'new description'

        The original biological sequence's ID and description have not been
        changed:

        >>> seq.id
        'id1'
        >>> seq.description
        'biological sequence'

        """
        defaults = {
            'sequence': self._bytes,
            'id': self.id,
            'description': self.description,
            'quality': self.quality
        }
        defaults.update(kwargs)
        return self._constructor(**defaults)

    def _constructor(self, **kwargs):
        return self.__class__(**kwargs)

    def equals(self, other, ignore=None, descriptive=False):
        """Compare two biological sequences for equality.

        By default, biological sequences are equal if their sequence,
        identifier, description, and quality scores are the same and they are
        the same type.

        Parameters
        ----------
        other : Sequence
            The sequence to test for equality against.
        ignore : iterable of str, optional
            List of features to ignore in the equality test. By default, all
            features must be the same for two biological sequences to be
            considered equal. Features that can be ignored are ``'type'``,
            ``'id'``, ``'description'``, ``'quality'``, and ``'sequence'``.

        Returns
        -------
        bool
            Indicates whether `self` and `other` are equal.

        See Also
        --------
        __eq__
        __ne__

        Examples
        --------
        Define two biological sequences that have the same underlying sequence
        of characters:

        >>> from skbio import Sequence
        >>> s = Sequence('GGUCGUGAAGGA')
        >>> t = Sequence('GGUCGUGAAGGA')

        The two sequences are considered equal because they are the same type,
        their underlying sequence of characters are the same, and their
        optional attributes (id, description, and quality scores) were not
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
        ``u``, but with different identifier and quality scores:
        >>> v = Sequence('GGUCGUGACCGA', id='abc',
        ...                        quality=[1, 5, 3, 3, 2, 42, 100, 9, 10, 55,
        ...                                 42, 42])

        By default, the two sequences are *not* considered equal because their
        identifiers and quality scores do not match:

        >>> u.equals(v)
        False

        By specifying that the quality scores and identifier should be ignored,
        they now compare equal:

        >>> u.equals(v, ignore=['quality', 'id'])
        True

        """
        if ignore is None:
            ignore = {}

        def raiser(feature, attr):
            if descriptive:
                raise AssertionError(
                    "%r is not equal to %r because of feature %r (%r != %r)" %
                    (self, other, feature,
                     getattr(self, attr), getattr(other, attr)))

        # Checks are ordered from least to most expensive.
        if 'type' not in ignore and self.__class__ != other.__class__:
            raiser('type', '__class__')
            return False

        if 'id' not in ignore and self.id != other.id:
            raiser('id', 'id')
            return False

        if 'description' not in ignore and \
                self.description != other.description:
            raiser('description', 'description')
            return False

        # Use array_equal instead of (a == b).all() because of this issue:
        #     http://stackoverflow.com/a/10582030
        if 'quality' not in ignore and not np.array_equal(self.quality,
                                                          other.quality):
            raiser('quality', 'quality')
            return False

        if 'sequence' not in ignore and self._string != other._string:
            raiser('sequence', 'sequence')
            return False

        return True

    def count(self, subsequence):
        """Returns the number of occurences of subsequence.

        Parameters
        ----------
        subsequence : str
            The subsequence to count occurences of.

        Returns
        -------
        int
            The number of occurrences of substring in the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUC')
        >>> s.count('G')
        2

        """
        if isinstance(subsequence, string_types):
            return self._string.count(subsequence)
        return self._string.count(Sequence(subsequence)._string)

    def distance(self, other, metric=None):
        """Returns the distance to other

        Parameters
        ----------
        other : `Sequence`
            The `Sequence` to compute the distance to.
        distance_fn : function, optional
            Function used to compute the distance between `self` and `other`.
            If ``None`` (the default), `scipy.spatial.distance.hamming` will be
            used.

        Returns
        -------
        float
            The distance between `self` and `other`.

        Raises
        ------
        skbio.sequence.SequenceError
            If ``len(self) != len(other)``

        See Also
        --------
        fraction_diff
        fraction_same
        skbio.DistanceMatrix
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('AGUC')
        >>> s.distance(t)
        0.25
        >>> def dumb_dist(s1, s2): return 0.42
        >>> s.distance(t, dumb_dist)
        0.42

        """
        # TODO refactor this method to accept a name (string) of the distance
        # metric to apply and accept **kwargs
        other = Sequence(other)
        if metric is None:
            # Hamming requires equal length sequences. We are checking this
            # here because the error you would get otherwise is cryptic.
            if len(self) != len(other):
                raise ValueError(
                    "Sequences do not have equal length. "
                    "Hamming distances can only be computed between "
                    "sequences of equal length.")
            metric = hamming
        return float(metric(self.sequence, other.sequence))

    def matches(self, other):
        if len(self) != len(other):
            raise ValueError("Match and mismatch vectors can only be "
                             "generated from equal length sequences.")
        return self._bytes == other._bytes

    def mismatches(self, other):
        return np.invert(self.matches(other))

    def mismatch_frequency(self, other, relative=False):
        """Return number of positions that differ relative to `other`

        Parameters
        ----------
        other : `Sequence`
            The `Sequence` to compare against.
        relative : ``bool``
            If ``True``, return the fraction of mismatches.

        Returns
        -------
        int, float
            The number of positions that differ between `self` and `other`.
            This will be an ``int`` if ``relative == False``, and a ``float``
            if ``relative == True``.

        Raises
        ------
        skbio.sequence.SequenceError
            If ``len(self) != len(other)``.

        See Also
        --------
        distance
        match_frequency
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('AGUC')
        >>> s.mismatch_frequency(t)
        1
        >>> s.mismatch_frequency(t, relative=True)
        0.25

        """
        if relative:
            return self.mismatches(other).mean()
        else:
            return self.mismatches(other).sum()

    def match_frequency(self, other, relative=False):
        """Return number of positions that are the same relative to `other`

        Parameters
        ----------
        other : `Sequence`
            The `Sequence` to compare against.
        relative : `bool`
            If ``True``, return the fraction of matches.

        Returns
        -------
        int, float
            The number of positions that are the same between `self` and
            `other`. This will be an ``int`` if ``relative == False``, and a
            ``float`` if ``relative == True``.

        Raises
        ------
        skbio.sequence.SequenceError
            If ``len(self) != len(other)``.

        See Also
        --------
        distance
        mismatch_frequency
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('GGUC')
        >>> t = Sequence('AGUC')
        >>> s.match_frequency(t)
        3
        >>> s.match_frequency(t, relative=True)
        0.75

        """
        if relative:
            return self.matches(other).mean()
        else:
            return self.matches(other).sum()

    def index(self, subsequence):
        """Return the position where subsequence first occurs

        Returns
        -------
        int
            The position where `subsequence` first occurs in the
            `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('ACACGACGTT-')
        >>> s.index('ACG')
        2

        """
        try:
            if isinstance(subsequence, string_types):
                return self._string.index(subsequence)
            return self._string.index(Sequence(subsequence)._string)
        except ValueError:
            raise ValueError(
                "%r is not present in %r." % (subsequence, self))

    def kmers(self, k, overlap=True):
        """Get the list of words of length k

        Parameters
        ----------
        k : int
            The word length.
        overlap : bool, optional
            Defines whether the k-words should be overlap or not
            overlap.

        Returns
        -------
        iterator of Sequences
            Iterator of words of length `k` contained in the
            Sequence.

        Raises
        ------
        ValueError
            If k < 1.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('ACACGACGTT')
        >>> [str(kw) for kw in s.k_words(4, overlap=False)]
        ['ACAC', 'GACG']
        >>> [str(kw) for kw in s.k_words(3, overlap=True)]
        ['ACA', 'CAC', 'ACG', 'CGA', 'GAC', 'ACG', 'CGT', 'GTT']

        """
        if k < 1:
            raise ValueError("k must be greater than 0.")

        if overlap:
            step = 1
        else:
            step = k

        for i in range(0, len(self) - k + 1, step):
            yield self[i:i+k]

    def kmer_frequencies(self, k, overlap=True, relative=False):
        """Get the counts of words of length k

        Parameters
        ----------
        k : int
            The word length.
        overlap : bool, optional
            Defines whether the k-words should be overlap or not
            overlap.

        Returns
        -------
        collections.Counter
            The counts of words of length `k` contained in the
            Sequence.

        Examples
        --------
        >>> from skbio.sequence import Sequence
        >>> s = Sequence('ACACAT')
        >>> s.k_word_counts(3, overlap=True)
        Counter({'ACA': 2, 'CAC': 1, 'CAT': 1})

        """
        kmers = self.kmers(k, overlap=overlap)
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

    def regex_iter(self, regex, retrieve_group_0=False):
        """Find patterns specified by regular expression

        Parameters
        ----------
        regex : SRE_Pattern
            A compiled regular expression (e.g., from re.compile) with
            finditer method
        retrieve_group_0 : bool, optional
            Defaults to ``False``. If ``True``, group(0) will be included in
            each list of tuples, which represents the shortest possible
            substring of the full sequence that contains all the other groups

        Returns
        -------
        generator
            yields lists of 3-tuples. Each 3-tuple represents a group from the
            matched regular expression, and contains the start of the hit, the
            end of the hit, and the substring that was hit
        """
        start = 0 if retrieve_group_0 else 1

        for match in regex.finditer(self._string):
            for g in range(start, len(match.groups())+1):
                yield slice(match.start(g), match.end(g))
