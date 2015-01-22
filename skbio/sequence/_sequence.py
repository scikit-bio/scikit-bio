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
from six import string_types

import re
import warnings
from collections import Sequence, Counter, defaultdict
from itertools import product

import numpy as np
from scipy.spatial.distance import hamming

from skbio._base import SkbioObject
from skbio.sequence import BiologicalSequenceError


class BiologicalSequence(Sequence, SkbioObject):
    """Base class for biological sequences.

    Parameters
    ----------
    sequence : python Sequence (e.g., str, list or tuple)
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
    validate : bool, optional
        If True, runs the `is_valid` method after construction and raises
        BiologicalSequenceError if ``is_valid == False``.

    Attributes
    ----------
    sequence
    id
    description
    quality

    Raises
    ------
    skbio.sequence.BiologicalSequenceError
        If ``validate == True`` and ``is_valid == False``, or if `quality` is
        not the correct shape.

    See Also
    --------
    NucleotideSequence
    DNASequence
    RNASequence

    Notes
    -----
    `BiologicalSequence` objects are immutable. Where applicable, methods
    return a new object of the same class.
    Subclasses are typically defined by methods relevant to only a specific
    type of biological sequence, and by containing characters only contained in
    the IUPAC standard character set [1]_ for that molecule type.

    Examples
    --------
    >>> from skbio.sequence import BiologicalSequence
    >>> s = BiologicalSequence('GGUCGUGAAGGA')
    >>> t = BiologicalSequence('GGUCCUGAAGGU')

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    """
    default_write_format = 'fasta'

    @classmethod
    def alphabet(cls):
        """Return the set of characters allowed in a `BiologicalSequence`.

        Returns
        -------
        set
            Characters that are allowed in a valid `BiologicalSequence`.

        See Also
        --------
        is_valid
        gap_alphabet
        unsupported_characters
        has_unsupported_characters

        """
        return cls.iupac_characters()

    @classmethod
    def gap_alphabet(cls):
        """Return the set of characters defined as gaps.

        Returns
        -------
        set
            Characters defined as gaps in a `BiologicalSequence`

        See Also
        --------
        alphabet
        unsupported_characters
        has_unsupported_characters
        degap
        gap_maps
        gap_vector

        """
        return set('-.')

    @classmethod
    def iupac_degenerate_characters(cls):
        """Return the degenerate IUPAC characters.

        Returns
        -------
        set
            Degenerate IUPAC characters.

        """
        return set(cls.iupac_degeneracies())

    @classmethod
    def iupac_characters(cls):
        """Return the non-degenerate and degenerate characters.

        Returns
        -------
        set
            Non-degenerate and degenerate characters.

        """
        return (cls.iupac_standard_characters() |
                cls.iupac_degenerate_characters())

    @classmethod
    def iupac_standard_characters(cls):
        """Return the non-degenerate IUPAC characters.

        Returns
        -------
        set
            Non-degenerate IUPAC characters.

        """
        return set()

    @classmethod
    def iupac_degeneracies(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate character to the set of
            non-degenerate IUPAC characters it represents.

        """
        return {}

    def __init__(self, sequence, id="", description="", quality=None,
                 validate=False):
        if not isinstance(sequence, string_types):
            sequence = ''.join(sequence)
        self._sequence = sequence

        self._id = id
        self._description = description
        self._set_quality(quality)

        if validate and not self.is_valid():
            unsupported_chars = self.unsupported_characters()
            raise BiologicalSequenceError(
                "Sequence contains unsupported characters: %s"
                % (" ".join(unsupported_chars)))

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
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> 'GGU' in s
        True
        >>> 'CCC' in s
        False

        .. shownumpydoc

        """
        return other in self._sequence

    def __eq__(self, other):
        """The equality operator.

        Biological sequences are equal if their sequence is the same and they
        are the same type. Identifier, description, and quality scores
        **are ignored**.

        Parameters
        ----------
        other : `BiologicalSequence`
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
        See ``BiologicalSequence.equals`` for more fine-grained control of
        equality testing.

        This method is equivalent to
        ``self.equals(other, ignore=['id', 'description', 'quality'])``.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> t = BiologicalSequence('GGUCGUGAAGGA')
        >>> s == t
        True
        >>> u = BiologicalSequence('GGUCGUGACCGA')
        >>> u == t
        False

        Note that even though the quality scores do not match between ``u`` and
        ``v``, they are considered equal:

        >>> v = BiologicalSequence('GGUCGUGACCGA',
        ...                        quality=[1, 5, 3, 3, 2, 42, 100, 9, 10, 55,
        ...                                 42, 42])
        >>> u == v
        True

        .. shownumpydoc

        """
        return self.equals(other, ignore=['id', 'description', 'quality'])

    def __getitem__(self, i):
        """The indexing operator.

        Parameters
        ----------
        i : int, slice, or sequence of ints
            The position(s) to return from the `BiologicalSequence`. If `i` is
            a sequence of ints, these are assumed to be indices in the sequence
            to keep.

        Returns
        -------
        BiologicalSequence
            New biological sequence containing the character(s) at position(s)
            `i` in the current `BiologicalSequence`. If quality scores are
            present, the quality score at position(s) `i` will be included in
            the returned sequence. ID and description are also included.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')

        Obtain a single character from the biological sequence:

        >>> s[1]
        <BiologicalSequence: G (length: 1)>

        Obtain a slice:

        >>> s[7:]
        <BiologicalSequence: AAGGA (length: 5)>

        Obtain characters at the following indices:

        >>> s[[3, 4, 7, 0, 3]]
        <BiologicalSequence: CGAGC (length: 5)>

        .. shownumpydoc

        """
        # TODO update this method when #60 is resolved. we have to deal with
        # discrepancies in indexing rules between str and ndarray... hence the
        # ugly code
        try:
            try:
                seq = self.sequence[i]
                qual = self.quality[i] if self.has_quality() else None
            except TypeError:
                seq = [self.sequence[idx] for idx in i]

                if self.has_quality():
                    qual = [self.quality[idx] for idx in i]
                else:
                    qual = None
        except IndexError:
            raise IndexError(
                "Position %r is out of range for %r." % (i, self))

        return self.copy(sequence=seq, quality=qual)

    def __hash__(self):
        """The hash operator.

        Returns
        -------
        int
            The hash of the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> hash(s)
        -1080059835405276950

        .. shownumpydoc

        """
        return hash(self._sequence)

    def __iter__(self):
        """The iter operator.

        Returns
        -------
        iterator
            Position iterator for the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> for c in s: print(c)
        G
        G
        U
        C

        .. shownumpydoc

        """
        return iter(self._sequence)

    def __len__(self):
        """The len operator.

        Returns
        -------
        int
            The length of the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> len(s)
        4

        .. shownumpydoc

        """
        return len(self._sequence)

    def __ne__(self, other):
        """The inequality operator.

        Biological sequences are not equal if their sequence is different or
        they are not the same type. Identifier, description, and quality scores
        **are ignored**.

        Parameters
        ----------
        other : `BiologicalSequence`
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
        See ``BiologicalSequence.equals`` for more fine-grained control of
        equality testing.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> t = BiologicalSequence('GGUCGUGAAGGA')
        >>> s != t
        False
        >>> u = BiologicalSequence('GGUCGUGACCGA')
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
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> repr(s)
        '<BiologicalSequence: GGUCGUGAAG... (length: 12)>'
        >>> t = BiologicalSequence('ACGT')
        >>> repr(t)
        '<BiologicalSequence: ACGT (length: 4)>'
        >>> t
        <BiologicalSequence: ACGT (length: 4)>

        .. shownumpydoc

        """
        first_ten = self.sequence[:10]
        cn = self.__class__.__name__
        length = len(self)
        if length > 10:
            ellipses = "..."
        else:
            ellipses = ""
        return '<%s: %s%s (length: %d)>' % (cn, first_ten, ellipses, length)

    def __reversed__(self):
        """The reversed operator.

        Returns
        -------
        iterator
            Reverse position iterator for the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> for c in reversed(s): print(c)
        C
        U
        G
        G

        .. shownumpydoc

        """
        return reversed(self._sequence)

    def __str__(self):
        """The str operator

        Returns
        -------
        str
            String representation of the `BiologicalSequence`. This will be the
            full sequence, but will not contain information about the type,
            identifier, description, or quality scores.

        See Also
        --------
        to_fasta
        id
        description
        __repr__

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> str(s)
        'GGUC'
        >>> print(s)
        GGUC

        .. shownumpydoc

        """
        return self.sequence

    @property
    def sequence(self):
        """String containing underlying biological sequence characters.

        A string representing the characters of the biological sequence.

        Notes
        -----
        This property is not writeable.

        """
        return self._sequence

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

    def has_quality(self):
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

    def copy(self, **kwargs):
        """Return a copy of the current biological sequence.

        Returns a copy of the current biological sequence, optionally with
        updated attributes specified as keyword arguments.

        Parameters
        ----------
        kwargs : dict, optional
            Keyword arguments passed to the ``BiologicalSequence`` (or
            subclass) constructor. The returned copy will have its attributes
            updated based on the values in `kwargs`. If an attribute is
            missing, the copy will keep the same attribute as the current
            biological sequence. Valid attribute names are `'sequence'`,
            `'id'`, `'description'`, and `'quality'`. Default behavior is to
            return a copy of the current biological sequence without changing
            any attributes.

        Returns
        -------
        BiologicalSequence
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

        >>> from skbio import BiologicalSequence
        >>> seq = BiologicalSequence('AACCGGTT', id='id1',
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
            'sequence': self.sequence,
            'id': self.id,
            'description': self.description,
            'quality': self.quality
        }
        defaults.update(kwargs)
        return self.__class__(**defaults)

    def equals(self, other, ignore=None):
        """Compare two biological sequences for equality.

        By default, biological sequences are equal if their sequence,
        identifier, description, and quality scores are the same and they are
        the same type.

        Parameters
        ----------
        other : BiologicalSequence
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

        >>> from skbio import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> t = BiologicalSequence('GGUCGUGAAGGA')

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

        >>> u = BiologicalSequence('GGUCGUGACCGA')
        >>> u.equals(t)
        False

        Define a biological sequence with the same sequence of characters as
        ``u``, but with different identifier and quality scores:
        >>> v = BiologicalSequence('GGUCGUGACCGA', id='abc',
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

        # Checks are ordered from least to most expensive.
        if 'type' not in ignore and self.__class__ != other.__class__:
            return False

        if 'id' not in ignore and self.id != other.id:
            return False

        if 'description' not in ignore and \
                self.description != other.description:
            return False

        # Use array_equal instead of (a == b).all() because of this issue:
        #     http://stackoverflow.com/a/10582030
        if 'quality' not in ignore and not np.array_equal(self.quality,
                                                          other.quality):
            return False

        if 'sequence' not in ignore and self.sequence != other.sequence:
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
            The number of occurrences of substring in the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> s.count('G')
        2

        """
        return self._sequence.count(subsequence)

    def degap(self):
        """Returns a new `BiologicalSequence` with gap characters removed.

        Returns
        -------
        BiologicalSequence
            A new `BiologicalSequence` with all characters from
            `self.gap_alphabet` filtered from the sequence.

        Notes
        -----
        The type, id, and description of the result will be the
        same as `self`. If quality scores are present, they will be filtered in
        the same manner as the sequence and included in the resulting
        degapped biological sequence.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC-C--ACGTT-C.', quality=range(16))
        >>> t = s.degap()
        >>> t
        <BiologicalSequence: GGUCCACGTT... (length: 11)>
        >>> print(t)
        GGUCCACGTTC
        >>> t.quality
        array([ 0,  1,  2,  3,  5,  8,  9, 10, 11, 12, 14])

        """
        gaps = self.gap_alphabet()
        indices = [i for i, e in enumerate(self) if e not in gaps]
        return self[indices]

    def distance(self, other, distance_fn=None):
        """Returns the distance to other

        Parameters
        ----------
        other : `BiologicalSequence`
            The `BiologicalSequence` to compute the distance to.
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
        skbio.sequence.BiologicalSequenceError
            If ``len(self) != len(other)`` and ``distance_fn`` ==
            ``scipy.spatial.distance.hamming``.

        See Also
        --------
        fraction_diff
        fraction_same
        skbio.DistanceMatrix
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> t = BiologicalSequence('AGUC')
        >>> s.distance(t)
        0.25
        >>> def dumb_dist(s1, s2): return 0.42
        >>> s.distance(t, dumb_dist)
        0.42

        """
        if distance_fn is None:
            distance_fn = hamming
            if len(self) != len(other):
                raise BiologicalSequenceError(
                    "Hamming distance can only be computed between "
                    "BiologicalSequences of equal length.")
        return distance_fn(self, other)

    def fraction_diff(self, other):
        """Return fraction of positions that differ relative to `other`

        Parameters
        ----------
        other : `BiologicalSequence`
            The `BiologicalSequence` to compare against.

        Returns
        -------
        float
            The fraction of positions that differ between `self` and `other`.

        Raises
        ------
        skbio.sequence.BiologicalSequenceError
            If ``len(self) != len(other)``.

        See Also
        --------
        distance
        fraction_same
        scipy.spatial.distance.hamming

        Notes
        -----
        Computed as the Hamming distance between `self` and `other`. This is
        available in addition to `distance` in case the `distance` method is
        updated to use something other than ``scipy.spatial.distance.hamming``
        as the default distance metric. So, if you specifically want the
        fraction of positions that differ, you should use this function instead
        of `distance` to ensure backward compatibility.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> t = BiologicalSequence('AGUC')
        >>> s.fraction_diff(t)
        0.25

        """
        return self.distance(other, distance_fn=hamming)

    def fraction_same(self, other):
        """Return fraction of positions that are the same relative to `other`

        Parameters
        ----------
        other : `BiologicalSequence`
            The `BiologicalSequence` to compare against.

        Returns
        -------
        float
            The fraction of positions that are the same between `self` and
            `other`.

        Raises
        ------
        skbio.sequence.BiologicalSequenceError
            If ``len(self) != len(other)``.

        See Also
        --------
        distance
        fraction_diff
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> t = BiologicalSequence('AGUC')
        >>> s.fraction_same(t)
        0.75

        """
        return 1. - self.fraction_diff(other)

    def gap_maps(self):
        """Return tuples mapping b/w gapped and ungapped positions

        Returns
        -------
        tuple containing two lists
            The first list is the length of the ungapped sequence, and each
            entry is the position of that base in the gapped sequence. The
            second list is the length of the gapped sequence, and each entry is
            either None (if that position represents a gap) or the position of
            that base in the ungapped sequence.

        See Also
        --------
        gap_vector

        Notes
        -----
        Visual aid is useful here. Imagine we have
        ``BiologicalSequence('-ACCGA-TA-')``. The position numbers in the
        ungapped sequence and gapped sequence will be as follows::

              0123456
              ACCGATA
              |||||\\
             -ACCGA-TA-
             0123456789

        So, in the first list, position 0 maps to position 1, position 1
        maps to position 2, position 5 maps to position 7, ... And, in the
        second list, position 0 doesn't map to anything (so it's None),
        position 1 maps to position 0, ...

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('-ACCGA-TA-')
        >>> m = s.gap_maps()
        >>> m[0]
        [1, 2, 3, 4, 5, 7, 8]
        >>> m[1]
        [None, 0, 1, 2, 3, 4, None, 5, 6, None]

        """
        degapped_to_gapped = []
        gapped_to_degapped = []
        non_gap_count = 0
        for i, e in enumerate(self):
            if self.is_gap(e):
                gapped_to_degapped.append(None)
            else:
                gapped_to_degapped.append(non_gap_count)
                degapped_to_gapped.append(i)
                non_gap_count += 1
        return degapped_to_gapped, gapped_to_degapped

    def gap_vector(self):
        """Return list indicating positions containing gaps

        Returns
        -------
        list of booleans
            The list will be of length ``len(self)``, and a position will
            contain ``True`` if the character at that position in the
            `BiologicalSequence` is in `self.gap_alphabet`, and ``False``
            otherwise.

        See Also
        --------
        gap_maps

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('..ACG--TT-')
        >>> s.gap_vector()
        [True, True, False, False, False, True, True, False, False, True]

        """
        return [self.is_gap(c) for c in self._sequence]

    def unsupported_characters(self):
        """Return the set of unsupported characters in the `BiologicalSequence`

        Returns
        -------
        set
            Invalid characters in the `BiologicalSequence` (i.e., the
            characters that are present in the `BiologicalSequence` but which
            are not in `BiologicalSequence.alphabet` or
            `BiologicalSequence.gap_alphabet`.

        See Also
        --------
        is_valid
        alphabet
        gap_alphabet
        has_unsupported_characters

        """
        return set(self) - self.alphabet() - self.gap_alphabet()

    def has_unsupported_characters(self):
        """Return bool indicating presence/absence of unsupported characters

        Returns
        -------
        bool
            ``True`` if invalid characters are present in the
            `BiologicalSequence` (i.e., characters which are not in
            `BiologicalSequence.alphabet` or
            `BiologicalSequence.gap_alphabet`) and ``False`` otherwise.

        See Also
        --------
        is_valid
        alphabet
        gap_alphabet
        has_unsupported_characters

        """
        all_supported = self.alphabet() | self.gap_alphabet()
        for e in self:
            if e not in all_supported:
                return True
        return False

    def index(self, subsequence):
        """Return the position where subsequence first occurs

        Returns
        -------
        int
            The position where `subsequence` first occurs in the
            `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('ACACGACGTT-')
        >>> s.index('ACG')
        2

        """
        try:
            return self._sequence.index(subsequence)
        except ValueError:
            raise ValueError(
                "%s is not present in %r." % (subsequence, self))

    @classmethod
    def is_gap(cls, char):
        """Return True if `char` is in the `gap_alphabet` set

        Parameters
        ----------
        char : str
            The string to check for presence in the `BiologicalSequence`
            `gap_alphabet`.

        Returns
        -------
        bool
            Indicates whether `char` is in the `BiologicalSequence` attribute
            `gap_alphabet`.

        Notes
        -----
        This is a class method.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> BiologicalSequence.is_gap('.')
        True
        >>> BiologicalSequence.is_gap('P')
        False
        >>> s = BiologicalSequence('ACACGACGTT')
        >>> s.is_gap('-')
        True

        """
        return char in cls.gap_alphabet()

    def is_gapped(self):
        """Return True if char(s) in `gap_alphabet` are present

        Returns
        -------
        bool
            Indicates whether there are one or more occurences of any character
            in `self.gap_alphabet` in the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('ACACGACGTT')
        >>> s.is_gapped()
        False
        >>> t = BiologicalSequence('A.CAC--GACGTT')
        >>> t.is_gapped()
        True

        """
        for e in self:
            if self.is_gap(e):
                return True
        return False

    def is_valid(self):
        """Return True if the sequence is valid

        Returns
        -------
        bool
            ``True`` if `self` is valid, and ``False`` otherwise.

        Notes
        -----
        Validity is defined as not containing any characters outside of
        `self.alphabet` and `self.gap_alphabet`.

        """
        return not self.has_unsupported_characters()

    def k_words(self, k, overlapping=True):
        """Get the list of words of length k

        Parameters
        ----------
        k : int
            The word length.
        overlapping : bool, optional
            Defines whether the k-words should be overlapping or not
            overlapping.

        Returns
        -------
        iterator of BiologicalSequences
            Iterator of words of length `k` contained in the
            BiologicalSequence.

        Raises
        ------
        ValueError
            If k < 1.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('ACACGACGTT')
        >>> [str(kw) for kw in s.k_words(4, overlapping=False)]
        ['ACAC', 'GACG']
        >>> [str(kw) for kw in s.k_words(3, overlapping=True)]
        ['ACA', 'CAC', 'ACG', 'CGA', 'GAC', 'ACG', 'CGT', 'GTT']

        """
        if k < 1:
            raise ValueError("k must be greater than 0.")

        sequence_length = len(self)

        if overlapping:
            step = 1
        else:
            step = k

        for i in range(0, sequence_length - k + 1, step):
            yield self[i:i+k]

    def k_word_counts(self, k, overlapping=True):
        """Get the counts of words of length k

        Parameters
        ----------
        k : int
            The word length.
        overlapping : bool, optional
            Defines whether the k-words should be overlapping or not
            overlapping.

        Returns
        -------
        collections.Counter
            The counts of words of length `k` contained in the
            BiologicalSequence.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('ACACAT')
        >>> s.k_word_counts(3, overlapping=True)
        Counter({'ACA': 2, 'CAC': 1, 'CAT': 1})

        """
        k_words = self.k_words(k, overlapping)
        return Counter((str(seq) for seq in k_words))

    def k_word_frequencies(self, k, overlapping=True):
        """Get the frequencies of words of length `k`

        Parameters
        ----------
        k : int
            The word length.
        overlapping : bool, optional
            Defines whether the k-words should be overlapping or not
            overlapping. This is only relevant when `k` > 1.

        Returns
        -------
        collections.defaultdict
            The frequencies of words of length `k` contained in the
            ``BiologicalSequence``.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('ACACAT')
        >>> s.k_word_frequencies(3, overlapping=True)
        defaultdict(<type 'float'>, {'CAC': 0.25, 'ACA': 0.5, 'CAT': 0.25})

        """
        if overlapping:
            num_words = len(self) - k + 1
        else:
            num_words = len(self) // k

        result = defaultdict(float)
        k_word_counts = self.k_word_counts(k, overlapping=overlapping)
        for word, count in viewitems(k_word_counts):
            result[str(word)] = count / num_words
        return result

    def lower(self):
        """Convert the BiologicalSequence to lowercase

        Returns
        -------
        BiologicalSequence
            The `BiologicalSequence` with all characters converted to
            lowercase.

        """
        return self.copy(sequence=self.sequence.lower())

    def nondegenerates(self):
        """Yield all nondegenerate versions of the sequence.

        Returns
        -------
        generator
            Generator yielding all possible nondegenerate versions of the
            sequence. Each sequence will have the same type, id, description,
            and quality scores as `self`.

        Raises
        ------
        BiologicalSequenceError
            If the sequence contains an invalid character (a character that
            isn't an IUPAC character or a gap character).

        See Also
        --------
        iupac_degeneracies

        Notes
        -----
        There is no guaranteed ordering to the generated sequences.

        Examples
        --------
        >>> from skbio.sequence import NucleotideSequence
        >>> seq = NucleotideSequence('TRG')
        >>> seq_generator = seq.nondegenerates()
        >>> for s in sorted(seq_generator, key=str): print(s)
        TAG
        TGG

        """
        degen_chars = self.iupac_degeneracies()
        nonexpansion_chars = self.iupac_standard_characters().union(
            self.gap_alphabet())

        expansions = []
        for char in self:
            if char in nonexpansion_chars:
                expansions.append(char)
            else:
                # Use a try/except instead of explicitly checking for set
                # membership on the assumption that an exception is rarely
                # thrown.
                try:
                    expansions.append(degen_chars[char])
                except KeyError:
                    raise BiologicalSequenceError(
                        "Sequence contains an invalid character: %s" % char)

        result = product(*expansions)
        return (self.copy(sequence=nondegen_seq) for nondegen_seq in result)

    def to_fasta(self, field_delimiter=" ", terminal_character="\n"):
        """Return the sequence as a fasta-formatted string

        .. note:: Deprecated in scikit-bio 0.2.0-dev
           ``to_fasta`` will be removed in scikit-bio 0.3.0. It is replaced by
           ``write``, which is a more general method for serializing
           FASTA-formatted files. ``write`` supports multiple file formats by
           taking advantage of scikit-bio's I/O registry system. See
           :mod:`skbio.io` for more details.

        Parameters
        ----------
        field_delimiter : str, optional
            The character(s) to use on the header line between the
            `self.id` and `self.description`.

        terminal_character : str, optional
            The last character to be included in the result (if you don't want
            a trailing newline or other character in the result, you can pass
            ``terminal_character=""``).

        Returns
        -------
        str
            The `BiologicalSequence` as a fasta-formatted string.

        Examples
        --------
        >>> from skbio.sequence import BiologicalSequence
        >>> s = BiologicalSequence('ACACGACGTT')
        >>> print(s.to_fasta(terminal_character=""))
        >
        ACACGACGTT
        >>> t = BiologicalSequence('ACA',id='my-seq',description='h')
        >>> print(t.to_fasta(terminal_character=""))
        >my-seq h
        ACA

        """
        warnings.warn(
            "BiologicalSequence.to_fasta is deprecated and will be removed in "
            "scikit-bio 0.3.0. Please update your code to use "
            "BiologicalSequence.write.", DeprecationWarning)

        if self._description:
            header_line = '%s%s%s' % (self._id, field_delimiter,
                                      self._description)
        else:
            header_line = self._id

        return '>%s\n%s%s' % (
            header_line, self.sequence, terminal_character)

    def upper(self):
        """Convert the BiologicalSequence to uppercase

        Returns
        -------
        BiologicalSequence
            The `BiologicalSequence` with all characters converted to
            uppercase.

        """
        return self.copy(sequence=self.sequence.upper())

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
                raise BiologicalSequenceError(
                    "Phred quality scores must be 1-D.")
            if len(quality) != len(self):
                raise BiologicalSequenceError(
                    "Number of Phred quality scores (%d) must match the "
                    "number of characters in the biological sequence (%d)." %
                    (len(quality), len(self._sequence)))
            if (quality < 0).any():
                raise BiologicalSequenceError(
                    "Phred quality scores must be greater than or equal to "
                    "zero.")

        self._quality = quality

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

        for match in regex.finditer(self._sequence):
            for g in range(start, len(match.groups())+1):
                yield (match.start(g), match.end(g), match.group(g))


class NucleotideSequence(BiologicalSequence):
    """Base class for nucleotide sequences.

    A `NucleotideSequence` is a `BiologicalSequence` with additional methods
    that are only applicable for nucleotide sequences, and containing only
    characters used in the IUPAC DNA or RNA lexicon.

    See Also
    --------
    BiologicalSequence

    Notes
    -----
    All uppercase and lowercase IUPAC DNA/RNA characters are supported.

    """

    @classmethod
    def complement_map(cls):
        """Return the mapping of characters to their complements.

        Returns
        -------
        dict
            Mapping of characters to their complements.

        Notes
        -----
        Complements cannot be defined for a generic `NucleotideSequence`
        because the complement of 'A' is ambiguous.
        `NucleotideSequence.complement_map` will therefore be the empty dict.
        Thanks, nature...

        """
        return {}

    @classmethod
    def iupac_standard_characters(cls):
        """Return the non-degenerate IUPAC nucleotide characters.

        Returns
        -------
        set
            Non-degenerate IUPAC nucleotide characters.

        """
        return set("ACGTUacgtu")

    @classmethod
    def iupac_degeneracies(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate nucleotide character to the set of
            non-degenerate IUPAC nucleotide characters it represents.

        """
        degen_map = {
            "R": set("AG"), "Y": set("CTU"), "M": set("AC"), "K": set("TUG"),
            "W": set("ATU"), "S": set("GC"), "B": set("CGTU"),
            "D": set("AGTU"), "H": set("ACTU"), "V": set("ACG"),
            "N": set("ACGTU")
        }

        for degen_char in list(degen_map.keys()):
            nondegen_chars = degen_map[degen_char]
            degen_map[degen_char.lower()] = set(
                ''.join(nondegen_chars).lower())

        return degen_map

    def _complement(self, reverse=False):
        """Returns `NucleotideSequence` that is (reverse) complement of `self`.

        Parameters
        ----------
        reverse : bool, optional
            If ``True``, reverse `self` before complementing.

        Returns
        -------
        NucelotideSequence
            The (reverse) complement of `self`. Specific type will be the same
            as ``type(self)``.

        Raises
        ------
        skbio.sequence.BiologicalSequenceError
            If a character is present in the `NucleotideSequence` that is not
            in the complement map.

        Notes
        -----
        This private method centralizes the logic for `complement` and
        `reverse_complement`.

        """
        result = []
        complement_map = self.complement_map()
        seq_iterator = reversed(self) if reverse else self
        for base in seq_iterator:
            try:
                result.append(complement_map[base])
            except KeyError:
                raise BiologicalSequenceError(
                    "Don't know how to complement base %s. Is it in "
                    "%s.complement_map?" % (base, self.__class__.__name__))

        quality = self.quality
        if self.has_quality() and reverse:
            quality = self.quality[::-1]

        return self.copy(sequence=result, quality=quality)

    def complement(self):
        """Return the complement of the `NucleotideSequence`

        Returns
        -------
        NucelotideSequence
            The complement of `self`. Specific type will be the same as
            ``type(self)``.

        Raises
        ------
        skbio.sequence.BiologicalSequenceError
            If a character is present in the `NucleotideSequence` that is not
            in `self.complement_map`.

        See Also
        --------
        reverse_complement
        complement_map

        Notes
        -----
        The type, id, description, and quality scores of the result will be the
        same as `self`.

        """
        return self._complement()

    def is_reverse_complement(self, other):
        """Return True if `other` is the reverse complement of `self`

        Returns
        -------
        bool
            `True` if `other` is the reverse complement of `self` and `False`
            otherwise.

        Raises
        ------
        skbio.sequence.BiologicalSequenceError
            If a character is present in `other` that is not in the
            `self.complement_map`.

        See Also
        --------
        reverse_complement

        """
        return self == other.reverse_complement()

    def reverse_complement(self):
        """Return the reverse complement of the `NucleotideSequence`

        Returns
        -------
        NucelotideSequence
            The reverse complement of `self`. Specific type will be the same as
            ``type(self)``.

        Raises
        ------
        skbio.sequence.BiologicalSequenceError
            If a character is present in the `NucleotideSequence` that is not
            in `self.complement_map`.

        See Also
        --------
        complement
        complement_map
        is_reverse_complement

        Notes
        -----
        The type, id, and description of the result will be the same as `self`.
        If quality scores are present, they will be reversed and included in
        the resulting biological sequence.

        """
        return self._complement(reverse=True)
    rc = reverse_complement

    def find_features(self, feature_type, min_length=1, allow_gaps=False):
        """Search the sequence for features

        Parameters
        ----------
        feature_type : {'purine_run', 'pyrimidine_run'}
            The type of feature to find
        min_length : int, optional
            Defaults to 1. Only features at least as long as this will be
            returned
        allow_gaps : bool, optional
            Defaults to ``False``. If ``True``, then gaps will not be
            considered to disrupt a feature

        Returns
        -------
        generator
            Yields tuples of the start of the feature, the end of the feature,
            and the subsequence that composes the feature

        Examples
        --------
        >>> from skbio.sequence import NucleotideSequence
        >>> s = NucleotideSequence('G-AT.T')
        >>> list(s.find_features('purine_run'))
        [(0, 1, 'G'), (2, 3, 'A')]
        >>> list(s.find_features('purine_run', 2))
        []
        >>> list(s.find_features('purine_run', 2, allow_gaps=True))
        [(0, 3, 'G-A')]
        >>> list(s.find_features('pyrimidine_run', 2, allow_gaps=True))
        [(3, 6, 'T.T')]

        """
        gaps = re.escape(''.join(self.gap_alphabet()))
        acceptable = gaps if allow_gaps else ''

        if feature_type == 'purine_run':
            pat_str = '([AGag%s]{%d,})' % (acceptable, min_length)
        elif feature_type == 'pyrimidine_run':
            pat_str = '([CTUctu%s]{%d,})' % (acceptable, min_length)
        else:
            raise ValueError("Unknown feature type: %s" % feature_type)

        pat = re.compile(pat_str)

        for hits in self.regex_iter(pat):
            if allow_gaps:
                degapped = hits[2]
                for gap_char in self.gap_alphabet():
                    degapped = degapped.replace(gap_char, '')
                if len(degapped) >= min_length:
                    yield hits
            else:
                yield hits


class DNASequence(NucleotideSequence):
    """Base class for DNA sequences.

    A `DNASequence` is a `NucelotideSequence` that is restricted to only
    containing characters used in IUPAC DNA lexicon.

    See Also
    --------
    NucleotideSequence
    BiologicalSequence

    Notes
    -----
    All uppercase and lowercase IUPAC DNA characters are supported.

    """

    @classmethod
    def complement_map(cls):
        """Return the mapping of characters to their complements.

        The complement of a gap character is itself.

        Returns
        -------
        dict
            Mapping of characters to their complements.

        """
        comp_map = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y',
            'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N', 'a': 't', 't': 'a', 'g': 'c',
            'c': 'g', 'y': 'r', 'r': 'y', 's': 's', 'w': 'w', 'k': 'm',
            'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', 'n': 'n'
        }

        comp_map.update({c: c for c in cls.gap_alphabet()})
        return comp_map

    @classmethod
    def iupac_standard_characters(cls):
        """Return the non-degenerate IUPAC DNA characters.

        Returns
        -------
        set
            Non-degenerate IUPAC DNA characters.

        """
        return set("ACGTacgt")

    @classmethod
    def iupac_degeneracies(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate DNA character to the set of
            non-degenerate IUPAC DNA characters it represents.

        """
        degen_map = {
            "R": set("AG"), "Y": set("CT"), "M": set("AC"), "K": set("TG"),
            "W": set("AT"), "S": set("GC"), "B": set("CGT"), "D": set("AGT"),
            "H": set("ACT"), "V": set("ACG"), "N": set("ACGT")
        }

        for degen_char in list(degen_map.keys()):
            nondegen_chars = degen_map[degen_char]
            degen_map[degen_char.lower()] = set(
                ''.join(nondegen_chars).lower())

        return degen_map


# class is accessible with alternative name for convenience
DNA = DNASequence


class RNASequence(NucleotideSequence):
    """Base class for RNA sequences.

    An `RNASequence` is a `NucelotideSequence` that is restricted to only
    containing characters used in the IUPAC RNA lexicon.

    Notes
    -----
    All uppercase and lowercase IUPAC RNA characters are supported.

    """

    @classmethod
    def complement_map(cls):
        """Return the mapping of characters to their complements.

        The complement of a gap character is itself.

        Returns
        -------
        dict
            Mapping of characters to their complements.

        """
        comp_map = {
            'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y',
            'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N', 'a': 'u', 'u': 'a', 'g': 'c',
            'c': 'g', 'y': 'r', 'r': 'y', 's': 's', 'w': 'w', 'k': 'm',
            'm': 'k', 'b': 'v', 'd': 'h', 'h': 'd', 'v': 'b', 'n': 'n'
        }

        comp_map.update({c: c for c in cls.gap_alphabet()})
        return comp_map

    @classmethod
    def iupac_standard_characters(cls):
        """Return the non-degenerate IUPAC RNA characters.

        Returns
        -------
        set
            Non-degenerate IUPAC RNA characters.

        """
        return set("ACGUacgu")

    @classmethod
    def iupac_degeneracies(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate RNA character to the set of
            non-degenerate IUPAC RNA characters it represents.

        """
        degen_map = {
            "R": set("AG"), "Y": set("CU"), "M": set("AC"), "K": set("UG"),
            "W": set("AU"), "S": set("GC"), "B": set("CGU"), "D": set("AGU"),
            "H": set("ACU"), "V": set("ACG"), "N": set("ACGU")
        }

        for degen_char in list(degen_map.keys()):
            nondegen_chars = degen_map[degen_char]
            degen_map[degen_char.lower()] = set(
                ''.join(nondegen_chars).lower())

        return degen_map

# class is accessible with alternative name for convenience
RNA = RNASequence


class ProteinSequence(BiologicalSequence):
    """Base class for protein sequences.

    A `ProteinSequence` is a `BiologicalSequence` containing only characters
    used in the IUPAC protein lexicon.

    See Also
    --------
    BiologicalSequence

    Notes
    -----
    All uppercase and lowercase IUPAC protein characters are supported.

    """

    @classmethod
    def iupac_standard_characters(cls):
        """Return the non-degenerate IUPAC protein characters.

        Returns
        -------
        set
            Non-degenerate IUPAC protein characters.

        """
        return set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")

    @classmethod
    def iupac_degeneracies(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate protein character to the set of
            non-degenerate IUPAC protein characters it represents.

        """
        degen_map = {
            "B": set("DN"), "Z": set("EQ"),
            "X": set("ACDEFGHIKLMNPQRSTVWY")
        }

        degen_map_lower = {}
        for degen_char in degen_map:
            nondegen_chars = degen_map[degen_char]
            degen_map_lower[degen_char.lower()] = set(
                ''.join(nondegen_chars).lower())

        degen_map.update(degen_map_lower)

        return degen_map

# class is accessible with alternative name for convenience
Protein = ProteinSequence
