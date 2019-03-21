# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from abc import ABCMeta, abstractproperty
from itertools import product
import re

import numpy as np

from skbio.util._decorator import (classproperty, overrides, stable,
                                   deprecated, experimental)
from skbio.util._misc import MiniRegistry
from ._sequence import Sequence


class GrammaredSequenceMeta(ABCMeta, type):
    def __new__(mcs, name, bases, dct):
        cls = super(GrammaredSequenceMeta, mcs).__new__(mcs, name, bases, dct)

        concrete_gap_chars = \
            type(cls.gap_chars) is not abstractproperty
        concrete_degenerate_map = \
            type(cls.degenerate_map) is not abstractproperty
        concrete_definite_chars = \
            type(cls.definite_chars) is not abstractproperty
        concrete_default_gap_char = \
            type(cls.default_gap_char) is not abstractproperty
        # degenerate_chars is not abstract but it depends on degenerate_map
        # which is abstract.
        concrete_degenerate_chars = concrete_degenerate_map

        # Only perform metaclass checks if none of the attributes on the class
        # are abstract.
        # TODO: Rather than hard-coding a list of attributes to check, we can
        # probably check all the attributes on the class and make sure none of
        # them are abstract.
        if (concrete_gap_chars and concrete_degenerate_map and
                concrete_definite_chars and concrete_default_gap_char and
                concrete_degenerate_chars):

            if cls.default_gap_char not in cls.gap_chars:
                raise TypeError(
                    "default_gap_char must be in gap_chars for class %s" %
                    name)

            if len(cls.gap_chars & cls.degenerate_chars) > 0:
                raise TypeError(
                    "gap_chars and degenerate_chars must not share any "
                    "characters for class %s" % name)

            for key in cls.degenerate_map.keys():
                for definite_char in cls.degenerate_map[key]:
                    if definite_char not in cls.definite_chars:
                        raise TypeError(
                            "degenerate_map must expand only to "
                            "characters included in definite_chars "
                            "for class %s" % name)

            if len(cls.degenerate_chars & cls.definite_chars) > 0:
                raise TypeError(
                    "degenerate_chars and definite_chars must not "
                    "share any characters for class %s" % name)

            if len(cls.gap_chars & cls.definite_chars) > 0:
                raise TypeError(
                    "gap_chars and definite_chars must not share any "
                    "characters for class %s" % name)

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


class GrammaredSequence(Sequence, metaclass=GrammaredSequenceMeta):
    """Store sequence data conforming to a character set.

    This is an abstract base class (ABC) that cannot be instantiated.

    This class is intended to be inherited from to create grammared sequences
    with custom alphabets.

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
    ...     def definite_chars(cls):
    ...         return set("ABC")
    ...
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
    --------------------------
    Stats:
        length: 8
        has gaps: False
        has degenerates: False
        has definites: True
    --------------------------
    0 ABABACAC

    >>> seq = CustomSequence('XXXXXX')
    >>> seq
    CustomSequence
    -------------------------
    Stats:
        length: 6
        has gaps: False
        has degenerates: True
        has definites: False
    -------------------------
    0 XXXXXX

    """
    __validation_mask = None
    __degenerate_codes = None
    __definite_char_codes = None
    __gap_codes = None

    @classproperty
    def _validation_mask(cls):
        # TODO These masks could be defined (as literals) on each concrete
        # object. For now, memoize!
        if cls.__validation_mask is None:
            as_bytes = ''.join(cls.alphabet).encode('ascii')
            cls.__validation_mask = np.invert(np.bincount(
                np.frombuffer(as_bytes, dtype=np.uint8),
                minlength=cls._number_of_extended_ascii_codes).astype(bool))
        return cls.__validation_mask

    @classproperty
    def _degenerate_codes(cls):
        if cls.__degenerate_codes is None:
            degens = cls.degenerate_chars
            cls.__degenerate_codes = np.asarray([ord(d) for d in degens])
        return cls.__degenerate_codes

    @classproperty
    def _definite_char_codes(cls):
        if cls.__definite_char_codes is None:
            definite_chars = cls.definite_chars
            cls.__definite_char_codes = np.asarray(
                [ord(d) for d in definite_chars])
        return cls.__definite_char_codes

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

        This includes gap, definite, and degenerate characters.

        Returns
        -------
        set
            Valid characters.

        """
        return cls.degenerate_chars | cls.definite_chars | cls.gap_chars

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
        raise NotImplementedError

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
        raise NotImplementedError

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

    @classproperty
    @deprecated(as_of='0.5.0', until='0.6.0',
                reason='Renamed to definite_chars')
    def nondegenerate_chars(cls):
        """Return non-degenerate characters.

        Returns
        -------
        set
            Non-degenerate characters.

        """
        return cls.definite_chars

    @abstractproperty
    @classproperty
    @stable(as_of='0.5.0')
    def definite_chars(cls):
        """Return definite characters.

        Returns
        -------
        set
            Definite characters.

        """
        raise NotImplementedError

    @abstractproperty
    @classproperty
    @stable(as_of='0.4.0')
    def degenerate_map(cls):
        """Return mapping of degenerate to definite characters.

        Returns
        -------
        dict (set)
            Mapping of each degenerate character to the set of
            definite characters it represents.

        """
        raise NotImplementedError

    @property
    def _motifs(self):
        return _motifs

    @overrides(Sequence)
    def __init__(self, sequence, metadata=None, positional_metadata=None,
                 interval_metadata=None, lowercase=False, validate=True):
        super(GrammaredSequence, self).__init__(
            sequence, metadata, positional_metadata,
            interval_metadata, lowercase)

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
        definites
        has_definites

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
        definites
        has_definites

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

    @stable(as_of='0.5.0')
    def definites(self):
        """Find positions containing definite characters in the sequence.

        Returns
        -------
        1D np.ndarray (bool)
            Boolean vector where ``True`` indicates a definite character
            is present at that position in the biological sequence.

        See Also
        --------
        has_definites
        degenerates

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACWGN')
        >>> s.definites()
        array([ True,  True, False,  True, False], dtype=bool)

        """
        return np.in1d(self._bytes, self._definite_char_codes)

    @deprecated(as_of='0.5.0', until='0.6.0',
                reason='Renamed to definites')
    def nondegenerates(self):
        """Find positions containing non-degenerate characters in the sequence.

        Returns
        -------
        1D np.ndarray (bool)
            Boolean vector where ``True`` indicates a non-degenerate character
            is present at that position in the biological sequence.

        See Also
        --------
        has_definites
        degenerates

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('ACWGN')
        >>> s.nondegenerates()
        array([ True,  True, False,  True, False], dtype=bool)

        """
        return self.definites()

    @stable(as_of='0.5.0')
    def has_definites(self):
        """Determine if sequence contains one or more definite characters

        Returns
        -------
        bool
            Indicates whether there are one or more occurrences of
            definite characters in the biological sequence.

        See Also
        --------
        definites
        degenerates
        has_degenerates

        Examples
        --------
        >>> from skbio import DNA
        >>> s = DNA('NWNNNNNN')
        >>> s.has_definites()
        False
        >>> t = DNA('ANCACWWGACGTT')
        >>> t.has_definites()
        True

        """
        # TODO: cache results
        return bool(self.definites().any())

    @deprecated(as_of='0.5.0', until='0.6.0',
                reason='Renamed to has_definites')
    def has_nondegenerates(self):
        """Determine if sequence contains one or more non-degenerate characters

        Returns
        -------
        bool
            Indicates whether there are one or more occurrences of
            non-degenerate characters in the biological sequence.

        See Also
        --------
        definites
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
        return self.has_definites()

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
            has definites: True
            GC-content: 55.56%
        -----------------------------
        0 GGTCCATTC

        """
        return self[np.invert(self.gaps())]

    @stable(as_of='0.4.0')
    def expand_degenerates(self):
        """Yield all possible definite versions of the sequence.

        Yields
        ------
        GrammaredSequence
            Definite version of the sequence.

        See Also
        --------
        degenerate_map

        Notes
        -----
        There is no guaranteed ordering to the definite sequences that are
        yielded.

        Each definite sequence will have the same type, metadata, and
        positional metadata as the biological sequence.

        Examples
        --------
        >>> from skbio import DNA
        >>> seq = DNA('TRG')
        >>> seq_generator = seq.expand_degenerates()
        >>> for s in sorted(seq_generator, key=str):
        ...     s
        ...     print('')
        DNA
        --------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 33.33%
        --------------------------
        0 TAG
        <BLANKLINE>
        DNA
        --------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 66.67%
        --------------------------
        0 TGG
        <BLANKLINE>

        """
        degen_chars = self.degenerate_map
        nonexpansion_chars = self.definite_chars.union(self.gap_chars)

        expansions = []
        for char in self:
            char = str(char)
            if char in nonexpansion_chars:
                expansions.append(char)
            else:
                expansions.append(degen_chars[char])

        metadata = None
        if self.has_metadata():
            metadata = self.metadata

        positional_metadata = None
        if self.has_positional_metadata():
            positional_metadata = self.positional_metadata

        for definite_seq in product(*expansions):
            yield self._constructor(
                sequence=''.join(definite_seq),
                metadata=metadata,
                positional_metadata=positional_metadata,
                interval_metadata=self.interval_metadata)

    @stable(as_of='0.4.1')
    def to_regex(self, within_capture=False):
        """Return regular expression object that accounts for degenerate chars.

        Parameters
        ----------
        within_capture : bool
            If ``True``, format the regex pattern for the sequence into a
            single capture group. If ``False``, compile the regex pattern as-is
            with no capture groups.

        Returns
        -------
        regex
            Pre-compiled regular expression object (as from ``re.compile``)
            that matches all definite versions of this sequence, and nothing
            else.

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
        >>> regex = seq.to_regex(within_capture=True)
        >>> regex.match('TAG').groups(0)
        ('TAG',)

        """
        regex_parts = []
        for base in str(self):
            if base in self.degenerate_chars:
                regex_parts.append('[{0}]'.format(
                    ''.join(self.degenerate_map[base])))
            else:
                regex_parts.append(base)

        regex_string = ''.join(regex_parts)

        if within_capture:
            regex_string = '({})'.format(regex_string)

        return re.compile(regex_string)

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
        stats.append(('has definites', '%r' % self.has_definites()))
        return stats


_motifs = MiniRegistry()

# Leave this at the bottom
_motifs.interpolate(GrammaredSequence, "find_motifs")
