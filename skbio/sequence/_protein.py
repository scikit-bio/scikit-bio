# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.util._decorator import classproperty, overrides
from ._grammared_sequence import GrammaredSequence, _motifs as parent_motifs


class Protein(GrammaredSequence):
    r"""Store protein sequence data and optional associated metadata.

    Parameters
    ----------
    sequence : str, Sequence, or 1D np.ndarray (np.uint8 or '\|S1')
        Characters representing the protein sequence itself.
    metadata : dict, optional
        Arbitrary metadata which applies to the entire sequence.
    positional_metadata : Pandas DataFrame consumable, optional
        Arbitrary per-character metadata. For example, quality data from
        sequencing reads. Must be able to be passed directly to the Pandas
        DataFrame constructor.
    interval_metadata : IntervalMetadata
        Arbitrary interval metadata which applies to intervals within
        a sequence to store interval features (such as protein domains).
    lowercase : bool or str, optional
        If ``True``, lowercase sequence characters will be converted to
        uppercase characters in order to be valid IUPAC Protein characters. If
        ``False``, no characters will be converted. If a str, it will be
        treated as a key into the positional metadata of the object. All
        lowercase characters will be converted to uppercase, and a ``True``
        value will be stored in a boolean array in the positional metadata
        under the key.
    validate : bool, optional
        If ``True``, validation will be performed to ensure that all sequence
        characters are in the IUPAC protein character set. If ``False``,
        validation will not be performed. Turning off validation will improve
        runtime performance. If invalid characters are present, however, there
        is **no guarantee that operations performed on the resulting object
        will work or behave as expected.** Only turn off validation if you are
        certain that the sequence characters are valid. To store sequence data
        that is not IUPAC-compliant, use ``Sequence``.

    See Also
    --------
    GrammaredSequence

    Notes
    -----
    According to the IUPAC notation [1]_ , a protein sequence may contain the
    following 20 definite characters (canonical amino acids):

    +-----+---------+--------------+
    |Code |3-letter |Amino acid    |
    +=====+=========+==============+
    |``A``|Ala      |Alanine       |
    +-----+---------+--------------+
    |``C``|Cys      |Cysteine      |
    +-----+---------+--------------+
    |``D``|Asp      |Aspartic acid |
    +-----+---------+--------------+
    |``E``|Glu      |Glutamic acid |
    +-----+---------+--------------+
    |``F``|Phe      |Phenylalanine |
    +-----+---------+--------------+
    |``G``|Gly      |Glycine       |
    +-----+---------+--------------+
    |``H``|His      |Histidine     |
    +-----+---------+--------------+
    |``I``|Ile      |Isoleucine    |
    +-----+---------+--------------+
    |``K``|Lys      |Lysine        |
    +-----+---------+--------------+
    |``L``|Leu      |Leucine       |
    +-----+---------+--------------+
    |``M``|Met      |Methionine    |
    +-----+---------+--------------+
    |``N``|Asn      |Asparagine    |
    +-----+---------+--------------+
    |``P``|Pro      |Proline       |
    +-----+---------+--------------+
    |``Q``|Gln      |Glutamine     |
    +-----+---------+--------------+
    |``R``|Arg      |Arginine      |
    +-----+---------+--------------+
    |``S``|Ser      |Serine        |
    +-----+---------+--------------+
    |``T``|Thr      |Threonine     |
    +-----+---------+--------------+
    |``V``|Val      |Valine        |
    +-----+---------+--------------+
    |``W``|Trp      |Tryptophan    |
    +-----+---------+--------------+
    |``Y``|Tyr      |Tyrosine      |
    +-----+---------+--------------+

    And the following four degenerate characters, each of which representing
    two or more amino acids:

    +-----+---------+------------+
    |Code |3-letter |Amino acids |
    +=====+=========+============+
    |``B``|Asx      |D or N      |
    +-----+---------+------------+
    |``Z``|Glx      |E or Q      |
    +-----+---------+------------+
    |``J``|Xle      |I or L      |
    +-----+---------+------------+
    |``X``|Xaa      |All 20      |
    +-----+---------+------------+

    Plus one stop character: ``*`` (Ter), and two gap characters: ``-`` and ``.``.

    Characters other than the above 27 are not allowed. If you intend to use
    additional characters to represent non-canonical amino acids, such as ``U``
    (Sec, Selenocysteine) and ``O`` (Pyl, Pyrrolysine), you may create a custom
    alphabet using ``GrammaredSequence``. Directly modifying the alphabet of
    ``Protein`` may break functions that rely on the IUPAC alphabet.

    It should be noted that some functions do not support certain characters.
    For example, the BLOSUM and PAM substitution matrices do not support ``J``
    (Xle). In such circumstances, unsupported characters will be replaced with
    ``X`` to represent any of the canonical amino acids.

    References
    ----------
    .. [1] Nomenclature for incompletely specified bases in nucleic acid
       sequences: recommendations 1984.
       Nucleic Acids Res. May 10, 1985; 13(9): 3021-3030.
       A Cornish-Bowden

    Examples
    --------
    >>> from skbio import Protein
    >>> Protein('PAW')
    Protein
    --------------------------
    Stats:
        length: 3
        has gaps: False
        has degenerates: False
        has definites: True
        has stops: False
    --------------------------
    0 PAW

    Convert lowercase characters to uppercase:

    >>> Protein('paW', lowercase=True)
    Protein
    --------------------------
    Stats:
        length: 3
        has gaps: False
        has degenerates: False
        has definites: True
        has stops: False
    --------------------------
    0 PAW

    """

    __stop_codes = None

    @classproperty
    def _stop_codes(cls):
        if cls.__stop_codes is None:
            stops = cls.stop_chars
            cls.__stop_codes = np.asarray([ord(s) for s in stops])
        return cls.__stop_codes

    @classproperty
    @overrides(GrammaredSequence)
    def alphabet(cls):
        return super(Protein, cls).alphabet | cls.stop_chars

    @classproperty
    @overrides(GrammaredSequence)
    def definite_chars(cls):
        return set("ACDEFGHIKLMNOPQRSTUVWY")

    @classproperty
    @overrides(GrammaredSequence)
    def noncanonical_chars(cls):
        return set("OU")

    @classproperty
    @overrides(GrammaredSequence)
    def degenerate_map(cls):
        return {
            "B": set("DN"),
            "Z": set("EQ"),
            "J": set("IL"),
            "X": set("ACDEFGHIKLMNOPQRSTUVWY"),
        }

    @classproperty
    def stop_chars(cls):
        """Return characters representing translation stop codons.

        Returns
        -------
        set
            Characters representing translation stop codons.

        """
        return set("*")

    @classproperty
    @overrides(GrammaredSequence)
    def gap_chars(cls):
        return set("-.")

    @classproperty
    @overrides(GrammaredSequence)
    def default_gap_char(cls):
        return "-"

    @classproperty
    @overrides(GrammaredSequence)
    def wildcard_char(cls):
        return "X"

    @property
    def _motifs(self):
        return _motifs

    def stops(self):
        """Find positions containing stop characters in the protein sequence.

        Returns
        -------
        1D np.ndarray (bool)
            Boolean vector where ``True`` indicates a stop character is present
            at that position in the protein sequence.

        See Also
        --------
        has_stops

        Examples
        --------
        >>> from skbio import Protein
        >>> s = Protein('PAW')
        >>> s.stops()
        array([False, False, False], dtype=bool)
        >>> s = Protein('PAW*E*')
        >>> s.stops()
        array([False, False, False,  True, False,  True], dtype=bool)

        """
        return np.in1d(self._bytes, self._stop_codes)

    def has_stops(self):
        """Determine if the sequence contains one or more stop characters.

        Returns
        -------
        bool
            Indicates whether there are one or more occurrences of stop
            characters in the protein sequence.

        Examples
        --------
        >>> from skbio import Protein
        >>> s = Protein('PAW')
        >>> s.has_stops()
        False
        >>> s = Protein('PAW*E*')
        >>> s.has_stops()
        True

        """
        return bool(self.stops().any())

    @overrides(GrammaredSequence)
    def _repr_stats(self):
        """Define custom statistics to display in the sequence's repr."""
        stats = super(Protein, self)._repr_stats()
        stats.append(("has stops", "%r" % self.has_stops()))
        return stats


_motifs = parent_motifs.copy()


@_motifs("N-glycosylation")
def _motif_nitro_glycosylation(sequence, min_length, ignore):
    """Identify N-glycosylation runs."""
    return sequence.find_with_regex("(N[^PX][ST][^PX])", ignore=ignore)


# Leave this at the bottom
_motifs.interpolate(Protein, "find_motifs")
