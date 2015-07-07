# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import skbio
from skbio.util import classproperty, overrides
from skbio.util._decorator import stable
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
        GC-content: 42.86%
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
        GC-content: 42.86%
    -----------------------------
    0 ACCGAAU

    """

    @classproperty
    @stable(as_of="0.4.0")
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
    @stable(as_of="0.4.0")
    @overrides(IUPACSequence)
    def nondegenerate_chars(cls):
        return set("ACGU")

    @classproperty
    @stable(as_of="0.4.0")
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

    @stable(as_of="0.4.0")
    def translate(self, genetic_code=1, *args, **kwargs):
        """Translate RNA sequence into protein sequence.

        Parameters
        ----------
        genetic_code : int, GeneticCode, optional
            Genetic code to use in translation. If ``int``, used as a table ID
            to look up the corresponding NCBI genetic code.
        args : tuple
            Positional arguments accepted by ``GeneticCode.translate``.
        kwargs : dict
            Keyword arguments accepted by ``GeneticCode.translate``.

        Returns
        -------
        Protein
            Translated sequence.

        See Also
        --------
        GeneticCode.translate
        GeneticCode.from_ncbi
        translate_six_frames

        Notes
        -----
        RNA sequence's metadata are included in the translated protein
        sequence. Positional metadata are not included.

        Examples
        --------
        Translate RNA into protein using NCBI's standard genetic code (table ID
        1, the default genetic code in scikit-bio):

        >>> from skbio import RNA
        >>> rna = RNA('AUGCCACUUUAA')
        >>> rna.translate()
        Protein
        -----------------------------
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            has stops: True
        -----------------------------
        0 MPL*

        Translate the same RNA sequence using a different NCBI genetic code
        (table ID 3, the yeast mitochondrial code) and specify that translation
        must terminate at the first stop codon:

        >>> rna.translate(3, stop='require')
        Protein
        -----------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            has stops: False
        -----------------------------
        0 MPT

        """
        if not isinstance(genetic_code, skbio.GeneticCode):
            genetic_code = skbio.GeneticCode.from_ncbi(genetic_code)
        return genetic_code.translate(self, *args, **kwargs)

    @stable(as_of="0.4.0")
    def translate_six_frames(self, genetic_code=1, *args, **kwargs):
        """Translate RNA into protein using six possible reading frames.

        The six possible reading frames are:

        * 1 (forward)
        * 2 (forward)
        * 3 (forward)
        * -1 (reverse)
        * -2 (reverse)
        * -3 (reverse)

        Translated sequences are yielded in this order.

        Parameters
        ----------
        genetic_code : int, GeneticCode, optional
            Genetic code to use in translation. If ``int``, used as a table ID
            to look up the corresponding NCBI genetic code.
        args : tuple
            Positional arguments accepted by
            ``GeneticCode.translate_six_frames``.
        kwargs : dict
            Keyword arguments accepted by ``GeneticCode.translate_six_frames``.

        Yields
        ------
        Protein
            Translated sequence in the current reading frame.

        See Also
        --------
        GeneticCode.translate_six_frames
        GeneticCode.from_ncbi
        translate

        Notes
        -----
        This method is faster than (and equivalent to) performing six
        independent translations using, for example:

        ``(seq.translate(reading_frame=rf)
        for rf in GeneticCode.reading_frames)``

        RNA sequence's metadata are included in each translated protein
        sequence. Positional metadata are not included.

        Examples
        --------
        Translate RNA into protein using the six possible reading frames and
        NCBI's standard genetic code (table ID 1, the default genetic code in
        scikit-bio):

        >>> from skbio import RNA
        >>> rna = RNA('AUGCCACUUUAA')
        >>> for protein in rna.translate_six_frames():
        ...     protein
        ...     print('')
        Protein
        -----------------------------
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            has stops: True
        -----------------------------
        0 MPL*
        <BLANKLINE>
        Protein
        -----------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            has stops: False
        -----------------------------
        0 CHF
        <BLANKLINE>
        Protein
        -----------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            has stops: False
        -----------------------------
        0 ATL
        <BLANKLINE>
        Protein
        -----------------------------
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            has stops: False
        -----------------------------
        0 LKWH
        <BLANKLINE>
        Protein
        -----------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            has stops: True
        -----------------------------
        0 *SG
        <BLANKLINE>
        Protein
        -----------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has non-degenerates: True
            has stops: False
        -----------------------------
        0 KVA
        <BLANKLINE>

        """
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
