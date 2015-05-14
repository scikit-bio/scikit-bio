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
from ._iupac_sequence import IUPACSequence, _motifs as parent_motifs


class NucleotideSequence(with_metaclass(ABCMeta, IUPACSequence)):
    """Abstract base class for storing an IUPAC nucleotide sequence.

    This is an abstract base class (ABC) that cannot be instantiated.

    Attributes
    ----------
    id
    description
    sequence
    quality
    alphabet
    gap_chars
    nondegenerate_chars
    degenerate_chars
    degenerate_map
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
            If ``True``, return the reverse complement. If quality scores are
            present, they will be reversed.

        Returns
        -------
        NucleotideSequence
            The (reverse) complement of the nucleotide sequence. The type, ID,
            description, and quality scores of the result will be the same as
            the nucleotide sequence. If `reverse` is ``True``, quality scores
            will be reversed if they are present.

        See Also
        --------
        reverse_complement
        complement_map

        Examples
        --------
        >>> from skbio import DNA
        >>> DNA('TTCATT', id='s', quality=range(6)).complement()
        DNA('AAGTAA', length=6, id='s', quality=[0, 1, 2, 3, 4, 5])
        >>> DNA('TTCATT', id='s', quality=range(6)).complement(reverse=True)
        DNA('AATGAA', length=6, id='s', quality=[5, 4, 3, 2, 1, 0])

        """
        result = self._complement_lookup[self._bytes]
        quality = self.quality
        if reverse:
            result = result[::-1]
            if self._has_quality():
                quality = self.quality[::-1]

        return self._to(sequence=result, quality=quality)

    def reverse_complement(self):
        """Return the reverse complement of the nucleotide sequence.

        Returns
        -------
        NucleotideSequence
            The reverse complement of the nucleotide sequence. The type, ID,
            and description of the result will be the same as the nucleotide
            sequence. If quality scores are present, they will be reversed.

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
        >>> DNA('TTCATT', id='s', quality=range(6)).reverse_complement()
        DNA('AATGAA', length=6, id='s', quality=[5, 4, 3, 2, 1, 0])

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


_motifs = parent_motifs.copy()


@_motifs("purine-run")
def _motif_purine_run(sequence, min_length, ignore):
    """Identifies purine runs"""
    return sequence.slices_from_regex("([AGR]{%d,})" % min_length,
                                      ignore=ignore)


@_motifs("pyrimidine-run")
def _motif_pyrimidine_run(sequence, min_length, ignore):
    """Identifies pyrimidine runs"""
    return sequence.slices_from_regex("([CTUY]{%d,})" % min_length,
                                      ignore=ignore)

# Leave this at the bottom
_motifs.interpolate(NucleotideSequence, "find_motifs")
