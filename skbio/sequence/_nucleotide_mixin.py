# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from abc import ABCMeta, abstractproperty

import numpy as np

from skbio.util._decorator import classproperty, stable
from ._grammared_sequence import _motifs as parent_motifs


class NucleotideMixin(metaclass=ABCMeta):
    """Mixin for adding funtionality for working with sequences of nucleotides.

    This is an abstract base class (ABC) that cannot be instantiated.

    See Also
    --------
    DNA
    RNA

    """
    __complement_lookup = None
    __gc_codes = None

    @classproperty
    def _complement_lookup(cls):
        if cls.__complement_lookup is not None:
            return cls.__complement_lookup

        lookup = np.zeros(cls._number_of_extended_ascii_codes, dtype=np.uint8)
        for key, value in cls.complement_map.items():
            lookup[ord(key)] = ord(value)
        cls.__complement_lookup = lookup
        return lookup

    @classproperty
    def _gc_codes(cls):
        if cls.__gc_codes is None:
            gc_iupac_chars = 'GCS'
            cls.__gc_codes = np.asarray([ord(g) for g in gc_iupac_chars])
        return cls.__gc_codes

    @property
    def _motifs(self):
        return _motifs

    @abstractproperty
    @classproperty
    @stable(as_of='0.4.0')
    def complement_map(cls):
        """Return mapping of nucleotide characters to their complements.

        Returns
        -------
        dict
            Mapping of each character to its complement.

        Notes
        -----
        Complements cannot be defined for a generic nucleotide sequence because
        the complement of ``A`` is ambiguous. Thanks, nature...

        """
        raise NotImplementedError

    @stable(as_of='0.4.0')
    def complement(self, reverse=False):
        """Return the complement of the nucleotide sequence.

        Parameters
        ----------
        reverse : bool, optional
            If ``True``, return the reverse complement. If positional and/or
            interval metadata are present, they will be reversed.

        Returns
        -------
        NucleotideMixin
            The (reverse) complement of the nucleotide sequence. The type and
            metadata of the result will be the same as the nucleotide
            sequence. If `reverse` is ``True``, positional or interval metadata
            will be reversed if it is present.

        See Also
        --------
        reverse_complement
        complement_map

        Examples
        --------
        >>> from skbio import DNA
        >>> seq = DNA('TTCATT', positional_metadata={'quality':range(6)})
        >>> seq
        DNA
        -----------------------------
        Positional metadata:
            'quality': <dtype: int64>
        Stats:
            length: 6
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 16.67%
        -----------------------------
        0 TTCATT
        >>> seq.complement()
        DNA
        -----------------------------
        Positional metadata:
            'quality': <dtype: int64>
        Stats:
            length: 6
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 16.67%
        -----------------------------
        0 AAGTAA
        >>> rc = seq.complement(reverse=True)
        >>> rc
        DNA
        -----------------------------
        Positional metadata:
            'quality': <dtype: int64>
        Stats:
            length: 6
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 16.67%
        -----------------------------
        0 AATGAA
        >>> rc.positional_metadata['quality'].values
        array([5, 4, 3, 2, 1, 0])

        """
        result = self._complement_lookup[self._bytes]

        metadata = None
        if self.has_metadata():
            metadata = self.metadata

        positional_metadata = None
        if self.has_positional_metadata():
            positional_metadata = self.positional_metadata

        complement = self._constructor(
            sequence=result,
            metadata=metadata,
            positional_metadata=positional_metadata)

        if reverse:
            # this has to be before the interval metadata code,
            # because __gititem__ drops interval_metadata.
            complement = complement[::-1]

        if self.has_interval_metadata():
            complement.interval_metadata = self.interval_metadata
            if reverse:
                # TODO: this can be revised to match
                # positional_metadata when __getitem__
                # supports interval_metadata
                complement.interval_metadata._reverse()

        return complement

    @stable(as_of='0.4.0')
    def reverse_complement(self):
        """Return the reverse complement of the nucleotide sequence.

        Returns
        -------
        NucleotideMixin
            The reverse complement of the nucleotide sequence. The type and
            metadata of the result will be the same as the nucleotide
            sequence. If positional metadata is present, it will be reversed.

        See Also
        --------
        complement
        is_reverse_complement

        Notes
        -----
        This method is equivalent to ``self.complement(reverse=True)``.

        Examples
        --------
        >>> from skbio import DNA
        >>> seq = DNA('TTCATT',
        ...           positional_metadata={'quality':range(6)})
        >>> seq = seq.reverse_complement()
        >>> seq
        DNA
        -----------------------------
        Positional metadata:
            'quality': <dtype: int64>
        Stats:
            length: 6
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 16.67%
        -----------------------------
        0 AATGAA
        >>> seq.positional_metadata['quality'].values
        array([5, 4, 3, 2, 1, 0])


        """
        return self.complement(reverse=True)

    @stable(as_of='0.4.0')
    def is_reverse_complement(self, other):
        r"""Determine if a sequence is the reverse complement of this sequence.

        Parameters
        ----------
        other : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
            Sequence to compare to.

        Returns
        -------
        bool
            ``True`` if `other` is the reverse complement of the nucleotide
            sequence.

        Raises
        ------
        TypeError
            If `other` is a ``Sequence`` object with a different type than the
            nucleotide sequence.

        See Also
        --------
        reverse_complement

        Examples
        --------
        >>> from skbio import DNA
        >>> DNA('TTCATT').is_reverse_complement('AATGAA')
        True
        >>> DNA('TTCATT').is_reverse_complement('AATGTT')
        False
        >>> DNA('ACGT').is_reverse_complement('ACGT')
        True

        """
        other = self._munge_to_sequence(other, 'is_reverse_complement')

        # avoid computing the reverse complement if possible
        if len(self) != len(other):
            return False
        else:
            # we reverse complement ourselves because `other` is a `Sequence`
            # object at this point and we only care about comparing the
            # underlying sequence data
            return self.reverse_complement()._string == other._string

    @stable(as_of='0.4.0')
    def gc_content(self):
        """Calculate the relative frequency of G's and C's in the sequence.

        This includes G, C, and S characters. This is equivalent to calling
        ``gc_frequency(relative=True)``. Note that the sequence will be
        degapped before the operation, so gap characters will not be included
        when calculating the length of the sequence.

        Returns
        -------
        float
            Relative frequency of G's and C's in the sequence.

        See Also
        --------
        gc_frequency

        Examples
        --------
        >>> from skbio import DNA
        >>> DNA('ACGT').gc_content()
        0.5
        >>> DNA('ACGTACGT').gc_content()
        0.5
        >>> DNA('ACTTAGTT').gc_content()
        0.25
        >>> DNA('ACGT--..').gc_content()
        0.5
        >>> DNA('--..').gc_content()
        0

        `S` means `G` or `C`, so it counts:

        >>> DNA('ASST').gc_content()
        0.5

        Other degenerates don't count:

        >>> DNA('RYKMBDHVN').gc_content()
        0.0

        """
        return self.gc_frequency(relative=True)

    @stable(as_of='0.4.0')
    def gc_frequency(self, relative=False):
        """Calculate frequency of G's and C's in the sequence.

        This calculates the minimum GC frequency, which corresponds to IUPAC
        characters G, C, and S (which stands for G or C).

        Parameters
        ----------
        relative : bool, optional
            If False return the frequency of G, C, and S characters (ie the
            count). If True return the relative frequency, ie the proportion
            of G, C, and S characters in the sequence. In this case the
            sequence will also be degapped before the operation, so gap
            characters will not be included when calculating the length of the
            sequence.

        Returns
        -------
        int or float
            Either frequency (count) or relative frequency (proportion),
            depending on `relative`.

        See Also
        --------
        gc_content

        Examples
        --------
        >>> from skbio import DNA
        >>> DNA('ACGT').gc_frequency()
        2
        >>> DNA('ACGT').gc_frequency(relative=True)
        0.5
        >>> DNA('ACGT--..').gc_frequency(relative=True)
        0.5
        >>> DNA('--..').gc_frequency(relative=True)
        0

        `S` means `G` or `C`, so it counts:

        >>> DNA('ASST').gc_frequency()
        2

        Other degenerates don't count:

        >>> DNA('RYKMBDHVN').gc_frequency()
        0

        """

        counts = np.bincount(self._bytes,
                             minlength=self._number_of_extended_ascii_codes)
        gc = counts[self._gc_codes].sum()
        if relative:
            seq = self.degap()
            if len(seq) != 0:
                gc /= len(seq)
        return gc


_motifs = parent_motifs.copy()


@_motifs("purine-run")
def _motif_purine_run(sequence, min_length, ignore):
    """Identifies purine runs"""
    return sequence.find_with_regex("([AGR]{%d,})" % min_length,
                                    ignore=ignore)


@_motifs("pyrimidine-run")
def _motif_pyrimidine_run(sequence, min_length, ignore):
    """Identifies pyrimidine runs"""
    return sequence.find_with_regex("([CTUY]{%d,})" % min_length,
                                    ignore=ignore)
