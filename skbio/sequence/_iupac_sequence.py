# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import with_metaclass

from abc import ABCMeta
from itertools import product

import numpy as np

from skbio.util import classproperty, overrides, abstractproperty, sphinx_hack
from skbio.util._misc import MiniRegistry
from ._sequence import Sequence


class IUPACSequence(with_metaclass(ABCMeta, Sequence)):
    """Store biological sequence data and optional associated metadata.

    Attributes
    ----------
    id
    description
    sequence
    quality
    alphabet
    nondegenerate_chars
    gap_chars
    degenerate_chars
    degenerate_map
    complement_map

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
    the IUPAC standard character set for that molecule type.

    Examples
    --------
    >>> from skbio.sequence import Sequence
    >>> s = Sequence('GGUCGUGAAGGA')
    >>> t = Sequence('GGUCCUGAAGGU')

    """
    _number_of_extended_ascii_codes = 256
    _ascii_lowercase_boundary = 90
    __validation_mask = None
    __degenerate_codes = None
    __nondegenerate_codes = None
    __gap_codes = None

    @classproperty
    def _validation_mask(cls):
        # TODO These masks could be defined (as literals) on each concrete
        # object. For now, memoize!
        if cls.__validation_mask is None:
            cls.__validation_mask = np.invert(np.bincount(
                np.fromstring(''.join(cls.alphabet), dtype=np.uint8),
                minlength=cls._number_of_extended_ascii_codes).astype(bool))
        return cls.__validation_mask

    @classproperty
    def _degenerate_codes(cls):
        if cls.__degenerate_codes is None:
            degens = cls.degenerate_chars
            cls.__degenerate_codes = np.asarray([ord(d) for d in degens])
        return cls.__degenerate_codes

    @classproperty
    def _nondegenerate_codes(cls):
        if cls.__nondegenerate_codes is None:
            nondegens = cls.nondegenerate_chars
            cls.__nondegenerate_codes = np.asarray([ord(d) for d in nondegens])
        return cls.__nondegenerate_codes

    @classproperty
    def _gap_codes(cls):
        if cls.__gap_codes is None:
            gaps = cls.gap_chars
            cls.__gap_codes = np.asarray([ord(g) for g in gaps])
        return cls.__gap_codes

    #@sphinx_hack
    @classproperty
    def alphabet(cls):
        """Return the allowed characters.

        Returns
        -------
        set
            Characters that are allowed in a valid sequence.

        .. shownumpydoc

        """
        return cls.degenerate_chars | cls.nondegenerate_chars | cls.gap_chars

    #@sphinx_hack
    @classproperty
    def gap_chars(cls):
        """Return the characters defined as gaps.

        Returns
        -------
        set
            Characters defined as gaps.

        .. shownumpydoc

        """
        return set('-.')

    #@sphinx_hack
    @classproperty
    def degenerate_chars(cls):
        """Return the degenerate IUPAC characters.

        Returns
        -------
        set
            Degenerate IUPAC characters.

        .. shownumpydoc

        """
        return set(cls.degenerate_map)

    @abstractproperty
    @classproperty
    def nondegenerate_chars(cls):
        """Return the non-degenerate IUPAC characters.

        Returns
        -------
        set
            Non-degenerate IUPAC characters.

        .. shownumpydoc

        """
        return set()  # pragma: no cover

    @abstractproperty
    @classproperty
    def degenerate_map(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate character to the set of
            non-degenerate IUPAC characters it represents.

        .. shownumpydoc

        """
        return set()  # pragma: no cover

    @property
    def _motifs(self):
            return _motifs

    @overrides(Sequence)
    def __init__(self, sequence, id="", description="", quality=None,
                 validate=True, case_insensitive=False):
        super(IUPACSequence, self).__init__(
            sequence, id=id, description=description, quality=quality)

        if case_insensitive:
            self._convert_to_uppercase()

        if validate:
            self._validate()

    def _convert_to_uppercase(self):
        lowercase = self._bytes > self._ascii_lowercase_boundary
        if np.any(lowercase):
            with self._byte_ownership():
                # ASCII is built such that the difference between uppercase and
                # lowercase is the 6th bit.
                self._bytes[lowercase] ^= 32

    def _validate(self):
        # This is the fastest way that we have found to identify the
        # presence or absence of certain characters (numbers).
        # It works by multiplying a mask where the numbers which are
        # permitted have a zero at their index, and all others have a one.
        # The result is a vector which will propogate counts of invalid
        # numbers and remove counts of valid numbers, so that we need only
        # see if the array is empty to determine validity.
        invalid_characters = np.bincount(
            self._bytes, minlength=self._number_of_extended_ascii_codes
        ) * self._validation_mask
        if np.any(invalid_characters):
            bad = list(np.where(
                invalid_characters > 0)[0].astype(np.uint8).view('|S1'))
            raise ValueError(
                "Invalid character%s in sequence: %r. Valid IUPAC characters: "
                "%r" % ('s' if len(bad) > 1 else '',
                        [str(b.tostring().decode("ascii")) for b in bad] if
                        len(bad) > 1 else bad[0],
                        list(self.alphabet)))

    def gaps(self):
        """Return vector indicating positions containing gaps

        Returns
        -------
        np.ndarray of bool
            Boolean vector where `True` indicates a gap character is present
            at that position in the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import DNA
        >>> s = DNA('AC-G-')
        >>> s.gaps()
        array([False, False,  True, False,  True], dtype=bool)

        .. shownumpydoc

        """
        return np.in1d(self._bytes, self._gap_codes)

    def has_gaps(self):
        """Return True if sequence contains one or more gap characters

        Returns
        -------
        bool
            Indicates whether there are one or more occurences of any character
            in `self.gap_chars` in the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import DNA
        >>> s = DNA('ACACGACGTT')
        >>> s.has_gaps()
        False
        >>> t = DNA('A.CAC--GACGTT')
        >>> t.has_gaps()
        True

        .. shownumpydoc

        """
        # TODO use count, there aren't that many gap chars
        return bool(self.gaps().any())

    def degenerates(self):
        """Return vector indicating positions containing a degenerate character

        Returns
        -------
        np.ndarray of bool
            Boolean vector where `True` indicates a degenerate character is
            present at that position in the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import DNA
        >>> s = DNA('ACWGN')
        >>> s.degenerates()
        array([False, False,  True, False,  True], dtype=bool)

        .. shownumpydoc

        """
        return np.in1d(self._bytes, self._degenerate_codes)

    def has_degenerates(self):
        """Return True if sequence contains one or more degenerate characters

        Returns
        -------
        bool
            Indicates whether there are one or more occurences of any character
            in `self.degenerate_chars` in the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import DNA
        >>> s = DNA('ACAC-GACGTT')
        >>> s.has_degenerates()
        False
        >>> t = DNA('ANCACWWGACGTT')
        >>> t.has_degenerates()
        True

        .. shownumpydoc

        """
        # TODO use bincount!
        return bool(self.degenerates().any())

    def nondegenerates(self):
        """Return vector indicating positions containing a non-degenerate char

        Returns
        -------
        np.ndarray of bool
            Boolean vector where `True` indicates a non-degenerate character is
            present at that position in the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import DNA
        >>> s = DNA('ACWGN')
        >>> s.nondegenerates()
        array([ True,  True, False,  True, False], dtype=bool)

        .. shownumpydoc

        """
        return np.in1d(self._bytes, self._nondegenerate_codes)

    def has_nondegenerates(self):
        """Return True if sequence contains one or more non-degenerate chars

        Returns
        -------
        bool
            Indicates whether there are one or more occurences of any character
            in `self.nondegenerate_chars` in the `Sequence`.

        Examples
        --------
        >>> from skbio.sequence import DNA
        >>> s = DNA('NWNNNNNN')
        >>> s.has_nondegenerates()
        False
        >>> t = DNA('ANCACWWGACGTT')
        >>> t.has_nondegenerates()
        True

        .. shownumpydoc

        """
        return bool(self.nondegenerates().any())

    def degap(self):
        """Return a new sequence with gap characters removed.

        Returns
        -------
        IUPACSequence
            A new sequence with all gap characters removed.

        See Also
        --------
        gap_chars

        Notes
        -----
        The type, id, and description of the result will be the
        same as `self`. If quality scores are present, they will be filtered in
        the same manner as the sequence and included in the resulting
        degapped sequence.

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('GGTC-C--ATT-C.', quality=range(14))
        >>> t = s.degap()
        >>> t
        DNA('GGTCCATTC', length=9, quality=[0, 1, 2, 3, 5, 8, 9, 10, 12])

        .. shownumpydoc

        """
        return self[np.invert(self.gaps())]

    def expand_degenerates(self):
        """Yield all nondegenerate versions of the sequence.

        Returns
        -------
        generator
            Generator yielding all possible nondegenerate versions of the
            sequence. Each sequence will have the same type, id, description,
            and quality scores as `self`.

        Raises
        ------
        SequenceError
            If the sequence contains an invalid character (a character that
            isn't an IUPAC character or a gap character).

        See Also
        --------
        degenerate_map

        Notes
        -----
        There is no guaranteed ordering to the generated sequences.

        Examples
        --------
        >>> from skbio import DNA
        >>> seq = DNA('TRG')
        >>> seq_generator = seq.expand_degenerates()
        >>> for s in sorted(seq_generator, key=str): print(s)
        TAG
        TGG

        .. shownumpydoc

        """
        degen_chars = self.degenerate_map
        nonexpansion_chars = self.nondegenerate_chars.union(self.gap_chars)

        expansions = []
        for char in self:
            char = str(char)
            if char in nonexpansion_chars:
                expansions.append(char)
            else:
                expansions.append(degen_chars[char])

        result = product(*expansions)
        return (self._to(sequence=''.join(nondegen_seq)) for nondegen_seq in
                result)

    def find_motifs(self, motif_type, min_length=1, exclude=None):
        """Search the sequence for motifs

        Options for `motif_type`:

        Parameters
        ----------
        motif_type : str
            The type of motif to find.
        min_length : int, optional
            Defaults to 1. Only motifs at least as long as this will be
            returned.
        exclude : array of bool
            A boolean vector indicating positions to ignore when matching.

        Returns
        -------
        generator
            Yields tuples of the start of the feature, the end of the feature,
            and the subsequence that composes the feature

        Raises
        ------
        ValueError
            If specifying an unknown motif type.

        Example
        -------
        >>> from skbio import DNA
        >>> s = DNA('ACGGGGAGGCGGAG')
        >>> for e in s.find_motifs('purine-run', min_length=2): print(e, s[e])
        slice(2, 9, None) GGGGAGG
        slice(10, 14, None) GGAG

        # Gap characters can disrupt motifs
        >>> s = DNA('GG-GG')
        >>> for e in s.find_motifs('purine-run'): print(e)
        slice(0, 2, None)
        slice(3, 5, None)

        # Gaps can be ignored by passing the gap vector to excludes:
        >>> s = DNA('GG-GG')
        >>> for e in s.find_motifs('purine-run', exclude=s.gaps()): print(e)
        slice(0, 5, None)

        .. shownumpydoc

        """
        if motif_type not in self._motifs:
            raise ValueError("Not a known motif (%r) for this sequence (%s)." %
                             (motif_type, self.__class__.__name__))

        return self._motifs[motif_type](self, min_length, exclude)

    @overrides(Sequence)
    def _constructor(self, **kwargs):
        return self.__class__(validate=False, case_insensitive=False, **kwargs)


_motifs = MiniRegistry()

# Leave this at the bottom
_motifs.interpolate(IUPACSequence, "find_motifs")
