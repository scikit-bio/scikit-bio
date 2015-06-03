# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util import classproperty, overrides
from ._nucleotide_sequence import NucleotideSequence
from ._iupac_sequence import IUPACSequence


class DNA(NucleotideSequence):
    """Store DNA sequence data and optional associated metadata.

    Only characters in the IUPAC DNA character set [1]_ are supported.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
        Characters representing the DNA sequence itself.
    id : str, optional
        Sequence identifier (e.g., an accession number).
    description : str, optional
        Description or comment about the sequence (e.g., "green fluorescent
        protein").
    quality : 1D array_like (int), optional
        Phred quality scores stored as nonnegative integers, one per sequence
        character. If provided, must be the same length as the DNA sequence.
        Can be a 1D ``np.ndarray`` of integers or a structure that can be
        converted into this representation using ``np.asarray``. A copy will
        *not* be made if `quality` is already a 1D ``np.ndarray`` with an
        ``int`` ``dtype``. The array will be made read-only (i.e., its
        ``WRITEABLE`` flag will be set to ``False``).
    validate : bool, optional
        If ``True``, validation will be performed to ensure that all sequence
        characters are in the IUPAC DNA character set. If ``False``, validation
        will not be performed. Turning off validation will improve runtime
        performance. If invalid characters are present, however, there is
        **no guarantee that operations performed on the resulting object will
        work or behave as expected.** Only turn off validation if you are
        certain that the sequence characters are valid. To store sequence data
        that is not IUPAC-compliant, use ``Sequence``.
    case_insenstive : bool, optional
        If ``True``, lowercase sequence characters will be converted to
        uppercase characters in order to be valid IUPAC DNA characters.

    Attributes
    ----------
    id
    description
    values
    quality
    alphabet
    gap_chars
    nondegenerate_chars
    degenerate_chars
    degenerate_map
    complement_map

    See Also
    --------
    RNA

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    Examples
    --------
    >>> from skbio import DNA
    >>> s = DNA('ACCGAAT')
    >>> s
    DNA('ACCGAAT', length=7)

    Convert lowercase characters to uppercase:

    >>> s = DNA('AcCGaaT', case_insensitive=True)
    >>> s
    DNA('ACCGAAT', length=7)

    """

    @classproperty
    @overrides(NucleotideSequence)
    def complement_map(cls):
        comp_map = {
            'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y',
            'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }

        comp_map.update({c: c for c in cls.gap_chars})
        return comp_map

    @classproperty
    @overrides(IUPACSequence)
    def nondegenerate_chars(cls):
        return set("ACGT")

    @classproperty
    @overrides(IUPACSequence)
    def degenerate_map(cls):
        return {
            "R": set("AG"), "Y": set("CT"), "M": set("AC"), "K": set("TG"),
            "W": set("AT"), "S": set("GC"), "B": set("CGT"), "D": set("AGT"),
            "H": set("ACT"), "V": set("ACG"), "N": set("ACGT")
        }
