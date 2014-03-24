#!/usr/bin/env python
r"""
Biological sequences (:mod:`skbio.core.sequence`)
=================================================

.. currentmodule:: skbio.core.sequence

This module provides functionality for working with biological sequences,
including generic sequences, nucelotide sequences, DNA sequences, and RNA
sequences. Class methods and attributes are also available to obtain valid
character sets, complement maps for different sequence types, and for
obtaining degenerate character definitions.

Classes
-------

.. autosummary::
   :toctree: generated/

   BiologicalSequence
   NucleotideSequence
   DNASequence
   RNASequence

Examples
--------
>>> from skbio.core.sequence import DNASequence, RNASequence

New sequences are created with optional identifier and description fields.

>>> d1 = DNASequence('ACC--G-GGTA..')
>>> d1 = DNASequence('ACC--G-GGTA..',identifier="seq1")
>>> d1 = DNASequence('ACC--G-GGTA..',identifier="seq1",description="GFP")

New sequences can also be created from existing sequences, for example as their
reverse complement or degapped (i.e., unaligned) version.

>>> d2 = d1.degap()
>>> d1
<DNASequence: ACC--G-GGT... (length: 13)>
>>> d2
<DNASequence: ACCGGGTA (length: 8)>
>>> d3 = d2.reverse_complement()
>>> d3
<DNASequence: TACCCGGT (length: 8)>

It's also straight-forward to compute distances between sequences (optionally
using user-defined distance metrics, default is Hamming distance) for use in
sequence clustering, phylogenetic reconstruction, etc.

>>> d4 = DNASequence('GACCCGCT')
>>> d5 = DNASequence('GACCCCCT')
>>> d3.distance(d4)
0.25
>>> d3.distance(d5)
0.375

Class-level methods contain information about the molecule types.

>>> DNASequence.iupac_degeneracies()['B']
set(['C', 'T', 'G'])

>>> RNASequence.iupac_degeneracies()['B']
set(['C', 'U', 'G'])

>>> DNASequence.is_gap('-')
True

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from collections import Sequence
from functools import total_ordering
from itertools import izip, product

from scipy.spatial.distance import hamming

from skbio.core.exception import BiologicalSequenceError


@total_ordering
class BiologicalSequence(Sequence):
    """Base class for biological sequences.

    Attributes
    ----------
    description
    identifier

    Notes
    -----
    `BiologicalSequence` objects are immutable. Where applicable, methods
    return a new object of the same class.
    """

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
        return set()

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

    def __init__(self, sequence, identifier="", description="",
                 validate=False):
        """Initialize a `BiologicalSequence` object.

        Parameters
        ----------
        sequence : python Sequence (e.g., str, list or tuple)
            The biological sequence.
        identifier : str, optional
            The sequence identifier (e.g., an accession number).
        description : str, optional
            A description or comment about the sequence (e.g., "green
            fluorescent protein").
        validate : bool, optional
            If True, runs the `is_valid` method after construction and raises
            BiologicalSequenceError if ``is_valid == False``.

        Raises
        ------
        skbio.core.exception.BiologicalSequenceError
          If ``validate == True`` and ``is_valid == False``.

        See Also
        --------
        NucleotideSequence
        DNASequence
        RNASequence

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> t = BiologicalSequence('GGUCCUGAAGGU')

        """
        self._sequence = ''.join(sequence)
        self._identifier = identifier
        self._description = description

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
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> 'GGU' in s
        True
        >>> 'CCC' in s
        False

        """
        return other in self._sequence

    def __eq__(self, other):
        """The equality operator.

        Parameters
        ----------
        other : `BiologicalSequence`
            The sequence to test for equality against.

        Returns
        -------
        bool
            Indicates whether `self` and `other` are equal.

        Notes
        -----
        `BiologicalSequences` are equal if their sequence is the same and
        they are the same type.

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> t = BiologicalSequence('GGUCGUGAAGGA')
        >>> s == t
        True
        >>> u = BiologicalSequence('GGUCGUGACCGA')
        >>> u == t
        False

        """
        if self.__class__ != other.__class__:
            return False
        elif self._sequence != other._sequence:
            return False
        else:
            return True

    def __getitem__(self, i):
        """The indexing operator.

        Parameters
        ----------
        i : int
            The position to return from the `BiologicalSequence`.

        Returns
        -------
        str
            The character at position `i` in the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> s[1]
        <BiologicalSequence: G (length: 1)>

        """
        try:
            return self.__class__(self._sequence[i],
                                  self.identifier, self.description)
        except IndexError:
            raise IndexError(
                "Position %d is out of range for %r." % (i, self))

    def __hash__(self):
        """The hash operator.

        Returns
        -------
        int
            The hash of the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> hash(s)
        -1080059835405276950

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
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> for c in s: print c
        G
        G
        U
        C

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
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> len(s)
        4

        """
        return len(self._sequence)

    def __lt__(self, other):
        return self._sequence < other._sequence

    def __ne__(self, other):
        """The inequality operator.

        Parameters
        ----------
        other : `BiologicalSequence`
            The sequence to test for inequality against.

        Returns
        -------
        bool
            Indicates whether `self` and `other` are not equal.

        Notes
        -----
        `BiologicalSequences` are not equal if their sequence is different or
        they are not the same type.

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> t = BiologicalSequence('GGUCGUGAAGGA')
        >>> s != t
        False
        >>> u = BiologicalSequence('GGUCGUGACCGA')
        >>> u != t
        True

        """
        return not self.__eq__(other)

    def __repr__(self):
        """The repr method.

        Returns
        -------
        str
            Returns a string representation of the object.

        Notes
        -----
        String representation contains the class name, the first ten
        characters of the sequence followed by elipses (or the full sequence
        and no elipses, if the sequence is less than 11 characters long),
        followed by the sequence length.

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUCGUGAAGGA')
        >>> repr(s)
        '<BiologicalSequence: GGUCGUGAAG... (length: 12)>'
        >>> t = BiologicalSequence('ACGT')
        >>> repr(t)
        '<BiologicalSequence: ACGT (length: 4)>'
        >>> t
        <BiologicalSequence: ACGT (length: 4)>

        """
        first_ten = str(self)[:10]
        cn = self.__class__.__name__
        length = len(self)
        if length > 10:
            elipses = "..."
        else:
            elipses = ""
        return '<%s: %s%s (length: %d)>' % (cn, first_ten, elipses, length)

    def __reversed__(self):
        """The reversed operator.

        Returns
        -------
        iterator
            Reverse position iterator for the `BiologicalSequence`.

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> for c in reversed(s): print c
        C
        U
        G
        G

        """
        return reversed(self._sequence)

    def __str__(self):
        """The str operator

        Returns
        -------
        str
            String representation of the `BiologicalSequence`. This will be the
            full sequence, but will not contain information about the type, or
            `self.identifier` or `self.description`.

        See Also
        --------
        to_fasta
        identifier
        description
        __repr__

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> str(s)
        'GGUC'
        >>> print s
        GGUC

        """
        return ''.join(self._sequence)

    @property
    def description(self):
        """Return the description of the `BiologicalSequence`

        Returns
        -------
        str
            The description attribute of the `BiologicalSequence`

        """
        return self._description

    @property
    def identifier(self):
        """Return the identifier of the `BiologicalSequence`

        Returns
        -------
        str
            The identifier attribute of the `BiologicalSequence`

        """
        return self._identifier

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
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC')
        >>> s.count('G')
        2

        """
        return self._sequence.count(subsequence)

    def degap(self):
        """Returns a new `BiologicalSequence` with gaps characters removed.

        Returns
        -------
        BiologicalSequence
            A new `BiologicalSequence` with all characters from
            `self.gap_alphabet` filtered from the sequence.

        Notes
        -----
        The type, identifier, and description of the result will be the
        same as `self`.

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('GGUC-C--ACGTT-C.')
        >>> t = s.degap()
        >>> t
        <BiologicalSequence: GGUCCACGTT... (length: 11)>
        >>> print t
        GGUCCACGTTC

        """
        gaps = self.gap_alphabet()
        result = [e for e in self._sequence if e not in gaps]
        return self.__class__(result, identifier=self._identifier,
                              description=self._description)

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
        skbio.core.exception.BiologicalSequenceError
            If ``len(self) != len(other)``.

        See Also
        --------
        fraction_diff
        fraction_same
        skbio.core.distance.DistanceMatrix
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
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
                "Distance can only be computed between BiologicalSequences "
                "of equal length.")
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
        skbio.core.exception.BiologicalSequenceError
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
        >>> from skbio.core.sequence import BiologicalSequence
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
        skbio.core.exception.BiologicalSequenceError
            If ``len(self) != len(other)``.

        See Also
        --------
        distance
        fraction_diff
        scipy.spatial.distance.hamming

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
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
        >>> from skbio.core.sequence import BiologicalSequence
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
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('..ACG--TT-')
        >>> s.gap_vector()
        [True, True, False, False, False, True, True, False, False, True]

        """
        return map(self.is_gap, self._sequence)

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
            if not e in all_supported:
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
        >>> from skbio.core.sequence import BiologicalSequence
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
        >>> from skbio.core.sequence import BiologicalSequence
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
        >>> from skbio.core.sequence import BiologicalSequence
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

    def lower(self):
        """Convert the BiologicalSequence to lowercase

        Returns
        -------
        BiologicalSequence
            The `BiologicalSequence` with all characters converted to
            lowercase.

        """
        return self.__class__(self._sequence.lower(),
                              self.identifier, self.description)

    def to_fasta(self, field_delimiter=" ", terminal_character="\n"):
        """Return the sequence as a fasta-formatted string

        Parameters
        ----------
        field_delimiter : str, optional
            The character(s) to use on the header line between the
            `self.identifier` and `self.description`.

        terminal_character : str, optional
            The last character to be included in the result (if you don't want
            a trailing newline or other character in the result, you can pass
            ``terminal_character=""``).

        Returns
        -------
        str
            The `BiologicalSequence` as a fasta-formatted string.

        See Also
        --------
        __str__

        Examples
        --------
        >>> from skbio.core.sequence import BiologicalSequence
        >>> s = BiologicalSequence('ACACGACGTT')
        >>> print s.to_fasta(terminal_character="")
        >
        ACACGACGTT
        >>> t = BiologicalSequence('ACA',identifier='my-seq',description='h')
        >>> print t.to_fasta(terminal_character="")
        >my-seq h
        ACA

        """
        if self._description:
            header_line = '%s%s%s' % (self._identifier, field_delimiter,
                                      self._description)
        else:
            header_line = self._identifier

        return '>%s\n%s%s' % (
            header_line, str(self), terminal_character)

    def upper(self):
        """Convert the BiologicalSequence to uppercase

        Returns
        -------
        BiologicalSequence
            The `BiologicalSequence` with all characters converted to
            uppercase.

        """
        return self.__class__(self._sequence.upper(),
                              self.identifier, self.description)


class NucleotideSequence(BiologicalSequence):
    """Base class for nucleotide sequences.

    Attributes
    ----------
    description
    identifier

    Notes
    -----
    All uppercase and lowercase IUPAC DNA/RNA characters are supported.

    """

    @classmethod
    def alphabet(cls):
        """Return the set of characters allowed in a `NucleotideSequence`.

        Returns
        -------
        set
            Characters that are allowed in a valid `NucleotideSequence`.

        """
        return cls.iupac_characters()

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

    @classmethod
    def iupac_degenerate_characters(cls):
        """Return the degenerate IUPAC nucleotide characters.

        Returns
        -------
        set
            Degenerate IUPAC nucleotide characters.

        """
        return set(cls.iupac_degeneracies())

    @classmethod
    def iupac_characters(cls):
        """Return the non-degenerate and degenerate nucleotide characters.

        Returns
        -------
        set
            Non-degenerate and degenerate nucleotide characters.

        """
        return (cls.iupac_standard_characters() |
                cls.iupac_degenerate_characters())

    def _complement(self, seq_iterator):
        """Returns `NucleotideSequence` that is complement of `seq_iterator`

        Parameters
        ----------
        seq_iterator : iterator
            The `BiologicalSequence` to be complemented.

        Returns
        -------
        NucelotideSequence
            The complement of the sequence represented by `seq_iterator`.
            Specific type will be the same as ``type(self)``.

        Raises
        ------
        skbio.core.exception.BiologicalSequenceError
            If a character is present in the `NucleotideSequence` that is not
            in the complement map.

        Notes
        -----
        This private method centralizes the logic for `complement` and
        `reverse_complement` by taking the sequence as an iterator (so it can
        be passed the result of either `iter` or `reversed`).

        """
        result = []
        complement_map = self.complement_map()
        for base in seq_iterator:
            try:
                result.append(complement_map[base])
            except KeyError:
                raise BiologicalSequenceError(
                    "Don't know how to complement base %s. Is it in "
                    "%s.complement_map?" % (base, self.__class__.__name__))
        return self.__class__(result, self._identifier, self._description)

    def complement(self):
        """Return the complement of the `NucleotideSequence`

        Returns
        -------
        NucelotideSequence
            The complement of `self`. Specific type will be the same as
            ``type(self)``.

        Raises
        ------
        skbio.core.exception.BiologicalSequenceError
            If a character is present in the `NucleotideSequence` that is not
            in `self.complement_map`.

        See Also
        --------
        reverse_complement
        complement_map

        """
        return self._complement(self)

    def is_reverse_complement(self, other):
        """Return True if `other` is the reverse complement of `self`

        Returns
        -------
        bool
            `True` if `other` is the reverse complement of `self` and `False`
            otherwise.

        Raises
        ------
        skbio.core.exception.BiologicalSequenceError
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
        skbio.core.exception.BiologicalSequenceError
            If a character is present in the `NucleotideSequence` that is not
            in `self.complement_map`.

        See Also
        --------
        complement
        complement_map
        is_reverse_complement

        """
        return self._complement(reversed(self))
    rc = reverse_complement

    def nondegenerates(self):
        """Yield all non-degenerate versions of the sequence."""
        # Should this return an empty list (as it does currently) or should it
        # return a list with an empty NucleotideSequence, i.e.
        # NucleotideSequence('') ?
        if len(self) < 1:
            return

        standard_chars = self.iupac_standard_characters()
        degen_chars = self.iupac_degeneracies()

        expansions = []
        for char in self:
            if char in standard_chars:
                expansions.append(set(char))
            else:
                expansions.append(degen_chars[char])

        result = product(*expansions)
        for nondegen_seq in result:
            yield self.__class__(''.join(nondegen_seq), self.identifier,
                                 self.description)


class DNASequence(NucleotideSequence):
    """Base class for DNA sequences.

    Attributes
    ----------
    description
    identifier

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

    Attributes
    ----------
    description
    identifier

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
