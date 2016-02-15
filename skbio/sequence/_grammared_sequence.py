# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from abc import ABCMeta, abstractproperty
from itertools import product
import re

import numpy as np
from six import add_metaclass


from skbio.util._decorator import (classproperty, overrides, stable,
                                   experimental)
from skbio.util._misc import MiniRegistry
from ._sequence import Sequence


class GrammaredSequenceException(TypeError):
    pass


class GrammaredSequenceMeta(ABCMeta, type):
    def __new__(mcs, name, bases, dct):
        cls = super(GrammaredSequenceMeta, mcs).__new__(mcs, name, bases, dct)

        # Only perform validation on classes that aren't abstract.
        if (type(cls.default_gap_char) is not abstractproperty and
                type(cls.gap_chars) is not abstractproperty):
            if cls.default_gap_char not in cls.gap_chars:
                raise GrammaredSequenceException(
                    "default_gap_char must be in gap_chars for class %s" %
                    name)

            if len(cls.gap_chars & cls.degenerate_chars) > 0:
                raise GrammaredSequenceException(
                    "gap_chars and degenerate_chars must not share any "
                    "characters for class %s" % name)

            if (type(cls.degenerate_map) is not abstractproperty and
                    type(cls.nondegenerate_chars) is not abstractproperty):

                for key in cls.degenerate_map.keys():
                    for nondegenerate in cls.degenerate_map[key]:
                        if nondegenerate not in cls.nondegenerate_chars:
                            raise GrammaredSequenceException(
                                "degenerate_map must expand only to "
                                "characters included in nondegenerate_chars "
                                "for class %s" % name)

                if len(cls.gap_chars & cls.nondegenerate_chars) > 0:
                    raise GrammaredSequenceException(
                        "gap_chars and nondegenerate_chars must not share any "
                        "characters for class %s" % name)

                if len(cls.degenerate_chars & cls.nondegenerate_chars) > 0:
                    raise GrammaredSequenceException(
                        "degenerate_chars and nondegenerate_chars must not "
                        "share any characters for class %s" % name)

        return cls


# Adapted from http://stackoverflow.com/a/16056691/943814
# Note that inheriting from GrammaredSequenceMeta, rather than something
# more general, is intentional. Multiple inheritance with metaclasses can be
# tricky and is not handled automatically in Python. Since this class needs to
# inherit both from ABCMeta and GrammaredSequenceMeta, the only way we could
# find to make this work was to have GrammaredSequenceMeta inherit from ABCMeta
# and then inherit from GrammaredSequenceMeta here.
class DisableSubclassingMeta(GrammaredSequenceMeta):
    def __new__(mcs, name, bases, dct):
        for b in bases:
            if isinstance(b, DisableSubclassingMeta):
                raise TypeError("Subclassing disabled for class %s. To create"
                                " a custom sequence class, inherit directly"
                                " from skbio.sequence.%s" %
                                (b.__name__, GrammaredSequence.__name__))
        return super(DisableSubclassingMeta, mcs).__new__(mcs, name, bases,
                                                          dict(dct))


@add_metaclass(GrammaredSequenceMeta)
class GrammaredSequence(Sequence):
    """Store sequence data conforming to a character set.

    This is an abstract base class (ABC) that cannot be instantiated.

    This class is intended to be inherited from to create grammared sequences
    with custom alphabets.

    Attributes
    ----------
    values
    metadata
    positional_metadata
    alphabet
    gap_chars
    default_gap_char
    nondegenerate_chars
    degenerate_chars
    degenerate_map

    Raises
    ------
    ValueError
        If sequence characters are not in the character set [1]_.

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

    Note in the example below that properties either need to be static or
    use skbio's `classproperty` decorator.

    >>> from skbio.sequence import GrammaredSequence
    >>> from skbio.util import classproperty
    >>> class CustomSequence(GrammaredSequence):
    ...     @classproperty
    ...     def degenerate_map(cls):
    ...         return {"X": set("AB")}
    ...
    ...     @classproperty
    ...     def nondegenerate_chars(cls):
    ...         return set("ABC")
    ...
    ...     @classproperty
    ...     def default_gap_char(cls):
    ...         return '-'
    ...
    ...     @classproperty
    ...     def gap_chars(cls):
    ...         return set('-.')

    >>> seq = CustomSequence('ABABACAC')
    >>> seq
    CustomSequence
    -----------------------------
    Stats:
        length: 8
        has gaps: False
        has degenerates: False
        has non-degenerates: True
    -----------------------------
    0 ABABACAC

    >>> seq = CustomSequence('XXXXXX')
    >>> seq
    CustomSequence
    ------------------------------
    Stats:
        length: 6
        has gaps: False
        has degenerates: True
        has non-degenerates: False
    ------------------------------
    0 XXXXXX

    """
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
    @stable(as_of='0.4.0')
    def alphabet(cls):
        """Return valid characters.

        This includes gap, non-degenerate, and degenerate characters.

        Returns
        -------
        set
            Valid characters.

        """
        return cls.degenerate_chars | cls.nondegenerate_chars | cls.gap_chars

    @abstractproperty
    @classproperty
    @stable(as_of='0.4.0')
    def gap_chars(cls):
        """Return characters defined as gaps.

        Returns
        -------
        set
            Characters defined as gaps.

        """
        return set('-.')

    @abstractproperty
    @classproperty
    @experimental(as_of='0.4.1')
    def default_gap_char(cls):
        """Gap character to use when constructing a new gapped sequence.

        This character is used when it is necessary to represent gap characters
        in a new sequence. For example, a majority consensus sequence will use
        this character to represent gaps.

        Returns
        -------
        str
            Default gap character.

        """
        return set()  # pragma: no cover

    @classproperty
    @stable(as_of='0.4.0')
    def degenerate_chars(cls):
        """Return degenerate characters.

        Returns
        -------
        set
            Degenerate characters.

        """
        return set(cls.degenerate_map)

    @abstractproperty
    @classproperty
    @stable(as_of='0.4.0')
    def nondegenerate_chars(cls):
        """Return non-degenerate characters.

        Returns
        -------
        set
            Non-degenerate characters.

        """
        return set()  # pragma: no cover

    @abstractproperty
    @classproperty
    @stable(as_of='0.4.0')
    def degenerate_map(cls):
        """Return mapping of degenerate to non-degenerate characters.

        Returns
        -------
        dict (set)
            Mapping of each degenerate character to the set of
            non-degenerate characters it represents.

        """
        return set()  # pragma: no cover

    @property
    def _motifs(self):
            return _motifs

    @overrides(Sequence)
    def __init__(self, sequence, metadata=None, positional_metadata=None,
                 lowercase=False, validate=True):
        super(GrammaredSequence, self).__init__(
            sequence, metadata, positional_metadata, lowercase)

        if validate:
            self._validate()

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
                "Invalid character%s in sequence: %r. \n"
                "Valid characters: %r\n"
                "Note: Use `lowercase` if your sequence contains lowercase "
                "characters not in the sequence's alphabet."
                % ('s' if len(bad) > 1 else '',
                   [str(b.tostring().decode("ascii")) for b in bad] if
                   len(bad) > 1 else bad[0],
                   list(self.alphabet)))

    @stable(as_of='0.4.0')
    def gaps(self):
        """Find positions containing gaps in the biological sequence.

        Returns
        -------
        1D np.ndarray (bool)
            Boolean vector where ``True`` indicates a gap character is present
            at that position in the biological sequence.

        See Also
        --------
        has_gaps

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('AC-G-')
        >>> s.gaps()
        array([False, False,  True, False,  True], dtype=bool)

        """
        return np.in1d(self._bytes, self._gap_codes)

    @stable(as_of='0.4.0')
    def has_gaps(self):
        """Determine if the sequence contains one or more gap characters.

        Returns
        -------
        bool
            Indicates whether there are one or more occurrences of gap
            characters in the biological sequence.

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACACGACGTT')
        >>> s.has_gaps()
        False
        >>> t = DNA('A.CAC--GACGTT')
        >>> t.has_gaps()
        True

        """
        # TODO use count, there aren't that many gap chars
        # TODO: cache results
        return bool(self.gaps().any())

    @stable(as_of='0.4.0')
    def degenerates(self):
        """Find positions containing degenerate characters in the sequence.

        Returns
        -------
        1D np.ndarray (bool)
            Boolean vector where ``True`` indicates a degenerate character is
            present at that position in the biological sequence.

        See Also
        --------
        has_degenerates
        nondegenerates
        has_nondegenerates

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACWGN')
        >>> s.degenerates()
        array([False, False,  True, False,  True], dtype=bool)

        """
        return np.in1d(self._bytes, self._degenerate_codes)

    @stable(as_of='0.4.0')
    def has_degenerates(self):
        """Determine if sequence contains one or more degenerate characters.

        Returns
        -------
        bool
            Indicates whether there are one or more occurrences of degenerate
            characters in the biological sequence.

        See Also
        --------
        degenerates
        nondegenerates
        has_nondegenerates

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACAC-GACGTT')
        >>> s.has_degenerates()
        False
        >>> t = DNA('ANCACWWGACGTT')
        >>> t.has_degenerates()
        True

        """
        # TODO use bincount!
        # TODO: cache results
        return bool(self.degenerates().any())

    @stable(as_of='0.4.0')
    def nondegenerates(self):
        """Find positions containing non-degenerate characters in the sequence.

        Returns
        -------
        1D np.ndarray (bool)
            Boolean vector where ``True`` indicates a non-degenerate character
            is present at that position in the biological sequence.

        See Also
        --------
        has_nondegenerates
        degenerates
        has_nondegenerates

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACWGN')
        >>> s.nondegenerates()
        array([ True,  True, False,  True, False], dtype=bool)

        """
        return np.in1d(self._bytes, self._nondegenerate_codes)

    @stable(as_of='0.4.0')
    def has_nondegenerates(self):
        """Determine if sequence contains one or more non-degenerate characters

        Returns
        -------
        bool
            Indicates whether there are one or more occurrences of
            non-degenerate characters in the biological sequence.

        See Also
        --------
        nondegenerates
        degenerates
        has_degenerates

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('NWNNNNNN')
        >>> s.has_nondegenerates()
        False
        >>> t = DNA('ANCACWWGACGTT')
        >>> t.has_nondegenerates()
        True

        """
        # TODO: cache results
        return bool(self.nondegenerates().any())

    @stable(as_of='0.4.0')
    def degap(self):
        """Return a new sequence with gap characters removed.

        Returns
        -------
        GrammaredSequence
            A new sequence with all gap characters removed.

        See Also
        --------
        gap_chars

        Notes
        -----
        The type and metadata of the result will be the same as the
        biological sequence. If positional metadata is present, it will be
        filtered in the same manner as the sequence characters and included in
        the resulting degapped sequence.

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('GGTC-C--ATT-C.',
        ...         positional_metadata={'quality':range(14)})
        >>> s.degap()
        DNA
        -----------------------------
        Positional metadata:
            'quality': <dtype: int64>
        Stats:
            length: 9
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            GC-content: 55.56%
        -----------------------------
        0 GGTCCATTC

        """
        return self[np.invert(self.gaps())]

    @stable(as_of='0.4.0')
    def expand_degenerates(self):
        """Yield all possible non-degenerate versions of the sequence.

        Yields
        ------
        GrammaredSequence
            Non-degenerate version of the sequence.

        See Also
        --------
        degenerate_map

        Notes
        -----
        There is no guaranteed ordering to the non-degenerate sequences that
        are yielded.

        Each non-degenerate sequence will have the same type, metadata,
        and positional metadata as the biological sequence.

        Examples
        --------
        >>> from skbio import DNA
        >>> seq = DNA('TRG')
        >>> seq_generator = seq.expand_degenerates()
        >>> for s in sorted(seq_generator, key=str):
        ...     s
        ...     print('')
        DNA
        -----------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            GC-content: 33.33%
        -----------------------------
        0 TAG
        <BLANKLINE>
        DNA
        -----------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            GC-content: 66.67%
        -----------------------------
        0 TGG
        <BLANKLINE>

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

    @stable(as_of='0.4.1')
    def to_regex(self):
        """Return regular expression object that accounts for degenerate chars.

        Returns
        -------
        regex
            Pre-compiled regular expression object (as from ``re.compile``)
            that matches all non-degenerate versions of this sequence, and
            nothing else.

        Examples
        --------
        >>> from skbio import DNA
        >>> seq = DNA('TRG')
        >>> regex = seq.to_regex()
        >>> regex.match('TAG').string
        'TAG'
        >>> regex.match('TGG').string
        'TGG'
        >>> regex.match('TCG') is None
        True

        """
        regex_string = []
        for base in str(self):
            if base in self.degenerate_chars:
                regex_string.append('[{0}]'.format(
                    ''.join(self.degenerate_map[base])))
            else:
                regex_string.append(base)
        return re.compile(''.join(regex_string))

    @stable(as_of='0.4.0')
    def find_motifs(self, motif_type, min_length=1, ignore=None):
        """Search the biological sequence for motifs.

        Options for `motif_type`:

        Parameters
        ----------
        motif_type : str
            Type of motif to find.
        min_length : int, optional
            Only motifs at least as long as `min_length` will be returned.
        ignore : 1D array_like (bool), optional
            Boolean vector indicating positions to ignore when matching.

        Yields
        ------
        slice
            Location of the motif in the biological sequence.

        Raises
        ------
        ValueError
            If an unknown `motif_type` is specified.

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACGGGGAGGCGGAG')
        >>> for motif_slice in s.find_motifs('purine-run', min_length=2):
        ...     motif_slice
        ...     str(s[motif_slice])
        slice(2, 9, None)
        'GGGGAGG'
        slice(10, 14, None)
        'GGAG'

        Gap characters can disrupt motifs:

        >>> s = DNA('GG-GG')
        >>> for motif_slice in s.find_motifs('purine-run'):
        ...     motif_slice
        slice(0, 2, None)
        slice(3, 5, None)

        Gaps can be ignored by passing the gap boolean vector to `ignore`:

        >>> s = DNA('GG-GG')
        >>> for motif_slice in s.find_motifs('purine-run', ignore=s.gaps()):
        ...     motif_slice
        slice(0, 5, None)

        """
        if motif_type not in self._motifs:
            raise ValueError("Not a known motif (%r) for this sequence (%s)." %
                             (motif_type, self.__class__.__name__))

        return self._motifs[motif_type](self, min_length, ignore)

    @overrides(Sequence)
    def _constructor(self, **kwargs):
        return self.__class__(validate=False, lowercase=False, **kwargs)

    @overrides(Sequence)
    def _repr_stats(self):
        """Define custom statistics to display in the sequence's repr."""
        stats = super(GrammaredSequence, self)._repr_stats()
        stats.append(('has gaps', '%r' % self.has_gaps()))
        stats.append(('has degenerates', '%r' % self.has_degenerates()))
        stats.append(('has non-degenerates', '%r' % self.has_nondegenerates()))
        return stats


_motifs = MiniRegistry()

# Leave this at the bottom
_motifs.interpolate(GrammaredSequence, "find_motifs")
