# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import with_metaclass

from abc import ABCMeta, abstractproperty
from itertools import product

import numpy as np

from skbio.util import classproperty, overrides
from ._sequence import Sequence


class IUPACSequence(with_metaclass(ABCMeta, Sequence)):
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

    @classproperty
    def alphabet(cls):
        """Return the set of characters allowed in an `IUPACSequence`.

        Returns
        -------
        set
            Characters that are allowed in a valid `IUPACSequence`.

        """
        return cls.degenerate_chars | cls.nondegenerate_chars | cls.gap_chars

    @classproperty
    def gap_chars(cls):
        """Return the set of characters defined as gaps.

        Returns
        -------
        set
            Characters defined as gaps.

        """
        return set('-.')

    @classproperty
    def degenerate_chars(cls):
        """Return the degenerate IUPAC characters.

        Returns
        -------
        set
            Degenerate IUPAC characters.

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

        """
        pass

    @abstractproperty
    @classproperty
    def degenerate_map(cls):
        """Return the mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict of sets
            Mapping of IUPAC degenerate character to the set of
            non-degenerate IUPAC characters it represents.

        """
        pass

    @overrides(Sequence)
    def __init__(self, sequence, id="", description="", quality=None,
                 validate=True, case_insensitive=False):
        super(IUPACSequence, self).__init__(sequence, id, description, quality)
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
            raise ValueError("Invalid character%s in sequence: %r" %
                             ('s' if len(bad) > 1 else '',
                              bad if len(bad) > 1 else bad[0]))

    @overrides(Sequence)
    def _constructor(self, **kwargs):
        return self.__class__(validate=False, case_insensitive=False, **kwargs)

    def gaps(self):
        return np.in1d(self._bytes, self._gap_codes)

    def has_gaps(self):
        """Return True if char(s) in `gap_chars` are present

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

        """
        # TODO use count, there aren't that many gap chars
        return bool(self.gaps().any())

    def degenerates(self):
        return np.in1d(self._bytes, self._degenerate_codes)

    def has_degenerates(self):
        # TODO use bincount!
        return bool(self.degenerates().any())

    def nondegenerates(self):
        return np.in1d(self._bytes, self._nondegenerate_codes)

    def has_nondegenerates(self):
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
