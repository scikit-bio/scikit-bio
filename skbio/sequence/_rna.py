# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

import skbio
from skbio.util import classproperty, overrides
from ._nucleotide_mixin import NucleotideMixin, _motifs as _parent_motifs
from ._iupac_sequence import IUPACSequence


class RNA(IUPACSequence, NucleotideMixin):
    """Store RNA sequence data and optional associated metadata.

    Only characters in the IUPAC RNA character set [1]_ are supported.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
        Characters representing the RNA sequence itself.
    metadata : dict, optional
        Arbitrary metadata which applies to the entire sequence.
    positional_metadata : Pandas DataFrame consumable, optional
        Arbitrary per-character metadata. For example, quality data from
        sequencing reads. Must be able to be passed directly to the Pandas
        DataFrame constructor.
    validate : bool, optional
        If ``True``, validation will be performed to ensure that all sequence
        characters are in the IUPAC RNA character set. If ``False``, validation
        will not be performed. Turning off validation will improve runtime
        performance. If invalid characters are present, however, there is
        **no guarantee that operations performed on the resulting object will
        work or behave as expected.** Only turn off validation if you are
        certain that the sequence characters are valid. To store sequence data
        that is not IUPAC-compliant, use ``Sequence``.
    case_insenstive : bool, optional
        If ``True``, lowercase sequence characters will be converted to
        uppercase characters in order to be valid IUPAC RNA characters.

    Attributes
    ----------
    values
    metadata
    positional_metadata
    alphabet
    gap_chars
    nondegenerate_chars
    degenerate_chars
    degenerate_map
    complement_map

    See Also
    --------
    DNA

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    Examples
    --------
    >>> from skbio import RNA
    >>> RNA('ACCGAAU')
    RNA
    -----------------------------
    Stats:
        length: 7
        has gaps: False
        has degenerates: False
        has non-degenerates: True
    -----------------------------
    0 ACCGAAU

    Convert lowercase characters to uppercase:

    >>> RNA('AcCGaaU', lowercase=True)
    RNA
    -----------------------------
    Stats:
        length: 7
        has gaps: False
        has degenerates: False
        has non-degenerates: True
    -----------------------------
    0 ACCGAAU

    """

    @classproperty
    @overrides(NucleotideMixin)
    def complement_map(cls):
        comp_map = {
            'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y',
            'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H',
            'H': 'D', 'V': 'B', 'N': 'N'
        }

        comp_map.update({c: c for c in cls.gap_chars})
        return comp_map

    @classproperty
    @overrides(IUPACSequence)
    def nondegenerate_chars(cls):
        return set("ACGU")

    @classproperty
    @overrides(IUPACSequence)
    def degenerate_map(cls):
        return {
            "R": set("AG"), "Y": set("CU"), "M": set("AC"), "K": set("UG"),
            "W": set("AU"), "S": set("GC"), "B": set("CGU"), "D": set("AGU"),
            "H": set("ACU"), "V": set("ACG"), "N": set("ACGU")
        }

    @property
    def _motifs(self):
        return _motifs

    def translate(self, genetic_code=1, *args, **kwargs):
        """

        Parameters
        ----------
        genetic_code : int, GeneticCode, optional
            The genetic code to use, will lookup the GeneticCode by int id or
            will use the GeneticCode provided.

        Returns
        -------
        Protein
            Translated protein sequence.
        """
        if not isinstance(genetic_code, skbio.GeneticCode):
            genetic_code = skbio.GeneticCode.from_ncbi(genetic_code)
        return genetic_code.translate(self, *args, **kwargs)

    def translate_six_frames(self, genetic_code=1, *args, **kwargs):
        if not isinstance(genetic_code, skbio.GeneticCode):
            genetic_code = skbio.GeneticCode.from_ncbi(genetic_code)
        return genetic_code.translate_six_frames(self, *args, **kwargs)

    @overrides(IUPACSequence)
    def _repr_stats(self):
        """Define custom statistics to display in the sequence's repr."""
        stats = super(RNA, self)._repr_stats()
        stats.append(('GC-content', '{:.2%}'.format(self.gc_content())))
        return stats


_motifs = _parent_motifs.copy()

# Leave this at the bottom
_motifs.interpolate(RNA, "find_motifs")
