# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import skbio
from skbio.util._decorator import classproperty, overrides
from ._nucleotide_mixin import NucleotideMixin, _motifs as _parent_motifs
from ._grammared_sequence import GrammaredSequence


class DNA(GrammaredSequence, NucleotideMixin):
    r"""Store DNA sequence data and optional associated metadata.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
        Characters representing the DNA sequence itself.
    metadata : dict, optional
        Arbitrary metadata which applies to the entire sequence.
    positional_metadata : Pandas DataFrame consumable, optional
        Arbitrary per-character metadata. For example, quality data from
        sequencing reads. Must be able to be passed directly to the Pandas
        DataFrame constructor.
    interval_metadata : IntervalMetadata
        Arbitrary interval metadata which applies to intervals within
        a sequence to store interval features (such as genes on the
        DNA sequence).
    lowercase : bool or str, optional
        If ``True``, lowercase sequence characters will be converted to
        uppercase characters in order to be valid IUPAC DNA characters. If
        ``False``, no characters will be converted. If a str, it will be
        treated as a key into the positional metadata of the object. All
        lowercase characters will be converted to uppercase, and a ``True``
        value will be stored in a boolean array in the positional metadata
        under the key.
    validate : bool, optional
        If ``True``, validation will be performed to ensure that all sequence
        characters are in the IUPAC DNA character set. If ``False``, validation
        will not be performed. Turning off validation will improve runtime
        performance. If invalid characters are present, however, there is
        **no guarantee that operations performed on the resulting object will
        work or behave as expected.** Only turn off validation if you are
        certain that the sequence characters are valid. To store sequence data
        that is not IUPAC-compliant, use ``Sequence``.

    See Also
    --------
    RNA
    GrammaredSequence

    Notes
    -----
    According to the IUPAC DNA character set [1]_ , a DNA sequence may contain
    the following four definite characters (canonical nucleotides):

    +-----+-----------+
    |Code |Nucleobase |
    +=====+===========+
    |``A``|Adenine    |
    +-----+-----------+
    |``C``|Cytosine   |
    +-----+-----------+
    |``G``|Guanine    |
    +-----+-----------+
    |``T``|Thymine    |
    +-----+-----------+

    And the following 11 degenerate characters, each of which representing 2-4
    nucleotides:

    +-----+-------------+-----------+
    |Code |Nucleobases  |Meaning    |
    +=====+=============+===========+
    |``R``|A or G       |Purine     |
    +-----+-------------+-----------+
    |``Y``|C or T       |Pyrimidine |
    +-----+-------------+-----------+
    |``S``|G or C       |Strong     |
    +-----+-------------+-----------+
    |``W``|A or T       |Weak       |
    +-----+-------------+-----------+
    |``K``|G or T       |Keto       |
    +-----+-------------+-----------+
    |``M``|A or C       |Amino      |
    +-----+-------------+-----------+
    |``B``|C, G or T    |Not A      |
    +-----+-------------+-----------+
    |``D``|A, G or T    |Not C      |
    +-----+-------------+-----------+
    |``H``|A, C or T    |Not G      |
    +-----+-------------+-----------+
    |``V``|A, C or G    |Not T      |
    +-----+-------------+-----------+
    |``N``|A, C, G or T |Any        |
    +-----+-------------+-----------+

    Plus two gap characters: ``-`` and ``.``.

    Characters other than the above 17 are not allowed. If you intend to use
    additional characters to represent non-canonical nucleobases, such as ``I``
    (Inosine), you may create a custom alphabet using ``GrammaredSequence``.
    Directly modifying the alphabet of ``DNA`` may break methods that rely on
    the IUPAC alphabet.

    It should be noted that some functions do not support degenerate characters
    characters. In such cases, they will be replaced with `N` to represent any
    of the canonical nucleotides.

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    Examples
    --------
    >>> from skbio import DNA
    >>> DNA('ACCGAAT')
    DNA
    --------------------------
    Stats:
        length: 7
        has gaps: False
        has degenerates: False
        has definites: True
        GC-content: 42.86%
    --------------------------
    0 ACCGAAT

    Convert lowercase characters to uppercase:

    >>> DNA('AcCGaaT', lowercase=True)
    DNA
    --------------------------
    Stats:
        length: 7
        has gaps: False
        has degenerates: False
        has definites: True
        GC-content: 42.86%
    --------------------------
    0 ACCGAAT

    """

    @classproperty
    @overrides(NucleotideMixin)
    def complement_map(cls):
        comp_map = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "Y": "R",
            "R": "Y",
            "S": "S",
            "W": "W",
            "K": "M",
            "M": "K",
            "B": "V",
            "D": "H",
            "H": "D",
            "V": "B",
            "N": "N",
        }

        comp_map.update({c: c for c in cls.gap_chars})
        return comp_map

    @classproperty
    @overrides(GrammaredSequence)
    def definite_chars(cls):
        return set("ACGT")

    @classproperty
    @overrides(GrammaredSequence)
    def degenerate_map(cls):
        return {
            "R": set("AG"),
            "Y": set("CT"),
            "M": set("AC"),
            "K": set("TG"),
            "W": set("AT"),
            "S": set("GC"),
            "B": set("CGT"),
            "D": set("AGT"),
            "H": set("ACT"),
            "V": set("ACG"),
            "N": set("ACGT"),
        }

    @classproperty
    @overrides(GrammaredSequence)
    def default_gap_char(cls):
        return "-"

    @classproperty
    @overrides(GrammaredSequence)
    def gap_chars(cls):
        return set("-.")

    @classproperty
    @overrides(GrammaredSequence)
    def wildcard_char(cls):
        return "N"

    @property
    def _motifs(self):
        return _motifs

    def transcribe(self):
        """Transcribe DNA into RNA.

        DNA sequence is assumed to be the coding strand. Thymine (T) is
        replaced with uracil (U) in the transcribed sequence.

        Returns
        -------
        RNA
            Transcribed sequence.

        See Also
        --------
        translate
        translate_six_frames

        Notes
        -----
        DNA sequence's metadata, positional, and interval
        metadata are included in the transcribed RNA sequence.

        Examples
        --------
        Transcribe DNA into RNA:

        >>> from skbio import DNA
        >>> dna = DNA('TAACGTTA')
        >>> dna
        DNA
        --------------------------
        Stats:
            length: 8
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 25.00%
        --------------------------
        0 TAACGTTA
        >>> dna.transcribe()
        RNA
        --------------------------
        Stats:
            length: 8
            has gaps: False
            has degenerates: False
            has definites: True
            GC-content: 25.00%
        --------------------------
        0 UAACGUUA

        """
        seq = self._string.replace(b"T", b"U")

        metadata = None
        if self.has_metadata():
            metadata = self.metadata

        positional_metadata = None
        if self.has_positional_metadata():
            positional_metadata = self.positional_metadata

        interval_metadata = None
        if self.has_interval_metadata():
            interval_metadata = self.interval_metadata

        # turn off validation because `seq` is guaranteed to be valid
        return skbio.RNA(
            seq,
            metadata=metadata,
            positional_metadata=positional_metadata,
            interval_metadata=interval_metadata,
            validate=False,
        )

    def translate(self, *args, **kwargs):
        """Translate DNA sequence into protein sequence.

        DNA sequence is assumed to be the coding strand. DNA sequence is first
        transcribed into RNA and then translated into protein.

        Parameters
        ----------
        args : tuple
            Positional arguments accepted by ``RNA.translate``.
        kwargs : dict
            Keyword arguments accepted by ``RNA.translate``.

        Returns
        -------
        Protein
            Translated sequence.

        See Also
        --------
        RNA.reverse_transcribe
        RNA.translate
        translate_six_frames
        transcribe

        Notes
        -----
        DNA sequence's metadata are included in the translated protein
        sequence. Positional metadata are not included.

        Examples
        --------
        Translate DNA into protein using NCBI's standard genetic code (table ID
        1, the default genetic code in scikit-bio):

        >>> from skbio import DNA
        >>> dna = DNA('ATGCCACTTTAA')
        >>> dna.translate()
        Protein
        --------------------------
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: True
        --------------------------
        0 MPL*

        Translate the same DNA sequence using a different NCBI genetic code
        (table ID 3, the yeast mitochondrial code) and specify that translation
        must terminate at the first stop codon:

        >>> dna.translate(3, stop='require')
        Protein
        --------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: False
        --------------------------
        0 MPT

        """
        return self.transcribe().translate(*args, **kwargs)

    def translate_six_frames(self, *args, **kwargs):
        """Translate DNA into protein using six possible reading frames.

        DNA sequence is assumed to be the coding strand. DNA sequence is first
        transcribed into RNA and then translated into protein. The six possible
        reading frames are:

        * 1 (forward)
        * 2 (forward)
        * 3 (forward)
        * -1 (reverse)
        * -2 (reverse)
        * -3 (reverse)

        Translated sequences are yielded in this order.

        Parameters
        ----------
        args : tuple
            Positional arguments accepted by ``RNA.translate_six_frames``.
        kwargs : dict
            Keyword arguments accepted by ``RNA.translate_six_frames``.

        Yields
        ------
        Protein
            Translated sequence in the current reading frame.

        See Also
        --------
        RNA.translate_six_frames
        translate
        transcribe

        Notes
        -----
        This method is faster than (and equivalent to) performing six
        independent translations using, for example:

        ``(seq.translate(reading_frame=rf)
        for rf in GeneticCode.reading_frames)``

        DNA sequence's metadata are included in each translated protein
        sequence. Positional metadata are not included.

        Examples
        --------
        Translate DNA into protein using the six possible reading frames and
        NCBI's standard genetic code (table ID 1, the default genetic code in
        scikit-bio):

        >>> from skbio import DNA
        >>> dna = DNA('ATGCCACTTTAA')
        >>> for protein in dna.translate_six_frames():
        ...     protein
        ...     print('')
        Protein
        --------------------------
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: True
        --------------------------
        0 MPL*
        <BLANKLINE>
        Protein
        --------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: False
        --------------------------
        0 CHF
        <BLANKLINE>
        Protein
        --------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: False
        --------------------------
        0 ATL
        <BLANKLINE>
        Protein
        --------------------------
        Stats:
            length: 4
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: False
        --------------------------
        0 LKWH
        <BLANKLINE>
        Protein
        --------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: True
        --------------------------
        0 *SG
        <BLANKLINE>
        Protein
        --------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: False
        --------------------------
        0 KVA
        <BLANKLINE>

        """
        return self.transcribe().translate_six_frames(*args, **kwargs)

    @overrides(GrammaredSequence)
    def _repr_stats(self):
        """Define custom statistics to display in the sequence's repr."""
        stats = super(DNA, self)._repr_stats()
        stats.append(("GC-content", "{:.2%}".format(self.gc_content())))
        return stats


_motifs = _parent_motifs.copy()

# Leave this at the bottom
_motifs.interpolate(DNA, "find_motifs")
