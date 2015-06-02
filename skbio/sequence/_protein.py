# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util import classproperty, overrides
from ._iupac_sequence import IUPACSequence, _motifs as parent_motifs


class Protein(IUPACSequence):
    """Store protein sequence data and optional associated metadata.

    Only characters in the IUPAC protein character set [1]_ are supported.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
        Characters representing the protein sequence itself.
    id : str, optional
        Sequence identifier (e.g., an accession number).
    description : str, optional
        Description or comment about the sequence (e.g., "green fluorescent
        protein").
    quality : 1D array_like (int), optional
        Phred quality scores stored as nonnegative integers, one per sequence
        character. If provided, must be the same length as the protein
        sequence. Can be a 1D ``np.ndarray`` of integers or a structure that
        can be converted into this representation using ``np.asarray``. A copy
        will *not* be made if `quality` is already a 1D ``np.ndarray`` with an
        ``int`` ``dtype``. The array will be made read-only (i.e., its
        ``WRITEABLE`` flag will be set to ``False``).
    validate : bool, optional
        If ``True``, validation will be performed to ensure that all sequence
        characters are in the IUPAC protein character set. If ``False``,
        validation will not be performed. Turning off validation will improve
        runtime performance. If invalid characters are present, however, there
        is **no guarantee that operations performed on the resulting object
        will work or behave as expected.** Only turn off validation if you are
        certain that the sequence characters are valid. To store sequence data
        that is not IUPAC-compliant, use ``Sequence``.
    case_insenstive : bool, optional
        If ``True``, lowercase sequence characters will be converted to
        uppercase characters in order to be valid IUPAC protein characters.

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

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    Examples
    --------
    >>> from skbio import Protein
    >>> s = Protein('PAW')
    >>> s
    Protein('PAW', length=3, has_metadata=False, has_positional_metadata=False)

    Convert lowercase characters to uppercase:

    >>> s = Protein('paW', case_insensitive=True)
    >>> s
    Protein('PAW', length=3, has_metadata=False, has_positional_metadata=False)

    """

    @property
    def _motifs(self):
        return _motifs

    @classproperty
    @overrides(IUPACSequence)
    def nondegenerate_chars(cls):
        return set("ACDEFGHIKLMNPQRSTVWY*")

    @classproperty
    @overrides(IUPACSequence)
    def degenerate_map(cls):
        return {
            "B": set("DN"), "Z": set("EQ"),
            "X": set("ACDEFGHIKLMNPQRSTVWY")
        }


_motifs = parent_motifs.copy()


@_motifs("N-glycosylation")
def _motif_nitro_glycosylation(sequence, min_length, ignore):
    """Identifies N-glycosylation runs"""
    return sequence.slices_from_regex("(N[^PX][ST][^PX])", ignore=ignore)

# Leave this at the bottom
_motifs.interpolate(Protein, "find_motifs")
