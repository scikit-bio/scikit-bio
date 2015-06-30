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

import numpy as np

from skbio.util import classproperty
from ._iupac_sequence import _motifs as parent_motifs


class NucleotideMixin(with_metaclass(ABCMeta, object)):
    """Mixin for adding funtionality for working with sequences of nucleotides.

    This is an abstract base class (ABC) that cannot be instantiated.

    Attributes
    ----------
    complement_map

    See Also
    --------
    DNA
    RNA

    """
    __complement_lookup = None

    @classproperty
    def _complement_lookup(cls):
        if cls.__complement_lookup is not None:
            return cls.__complement_lookup

        lookup = np.zeros(cls._number_of_extended_ascii_codes, dtype=np.uint8)
        for key, value in cls.complement_map.items():
            lookup[ord(key)] = ord(value)
        cls.__complement_lookup = lookup
        return lookup

    @property
    def _motifs(self):
        return _motifs

    @abstractproperty
    @classproperty
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
        return set()  # pragma: no cover

    def complement(self, reverse=False):
        """Return the complement of the nucleotide sequence.

        Parameters
        ----------
        reverse : bool, optional
            If ``True``, return the reverse complement. If positional metadata
            is present, it will be reversed.

        Returns
        -------
        NucleotideMixin
            The (reverse) complement of the nucleotide sequence. The type and
            metadata of the result will be the same as the nucleotide
            sequence. If `reverse` is ``True``, positional metadata
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
            has non-degenerates: True
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
            has non-degenerates: True
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
            has non-degenerates: True
        -----------------------------
        0 AATGAA
        >>> rc.positional_metadata['quality'].values
        array([5, 4, 3, 2, 1, 0])

        """
        result = self._complement_lookup[self._bytes]
        complement = self._to(sequence=result)
        if reverse:
            complement = complement[::-1]
        return complement

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
            has non-degenerates: True
        -----------------------------
        0 AATGAA
        >>> seq.positional_metadata['quality'].values
        array([5, 4, 3, 2, 1, 0])


        """
        return self.complement(reverse=True)

    def is_reverse_complement(self, other):
        """Determine if a sequence is the reverse complement of this sequence.

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

    def gc_content(self):
        """Calculate the relative frequency of G's and C's in the sequence.
        
        This includes G, C, and S characters. This is equivalent to calling
        ``gc_frequency(relative=True)``.

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

        """
        return self.gc_frequency(relative=True)

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

        """

        counts = np.bincount(self._bytes,
                             minlength=self._number_of_extended_ascii_codes)
        gc_ord = (ord(y) for y in 'CGS')
        gc = sum(counts[x] for x in gc_ord)
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
