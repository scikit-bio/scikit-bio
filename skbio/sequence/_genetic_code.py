# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.util._decorator import classproperty, classonlymethod
from skbio._base import SkbioObject
from skbio.sequence import Protein, RNA
from skbio._base import ElasticLines


class GeneticCode(SkbioObject):
    """Genetic code for translating codons to amino acids.

    Parameters
    ----------
    amino_acids : consumable by ``skbio.Protein`` constructor
        64-character vector containing IUPAC amino acid characters. The order
        of the amino acids should correspond to NCBI's codon order (see *Notes*
        section below). `amino_acids` is the "AAs" field in NCBI's genetic
        code format [1]_.
    starts : consumable by ``skbio.Protein`` constructor
        64-character vector containing only M and - characters, with start
        codons indicated by M. The order of the amino acids should correspond
        to NCBI's codon order (see *Notes* section below). `starts` is the
        "Starts" field in NCBI's genetic code format [1]_.
    name : str, optional
        Genetic code name. This is simply metadata and does not affect the
        functionality of the genetic code itself.

    See Also
    --------
    RNA.translate
    DNA.translate
    GeneticCode.from_ncbi

    Notes
    -----
    The genetic codes available via ``GeneticCode.from_ncbi`` and used
    throughout the examples are defined in [1]_. The genetic code strings
    defined there are directly compatible with the ``GeneticCode`` constructor.

    The order of `amino_acids` and `starts` should correspond to NCBI's codon
    order, defined in [1]_::

        UUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
        UUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGG
        UCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAG

    Note that scikit-bio displays this ordering using the IUPAC RNA alphabet,
    while NCBI displays this same ordering using the IUPAC DNA alphabet (for
    historical purposes).

    References
    ----------
    .. [1] http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

    Examples
    --------
    Get NCBI's standard genetic code (table ID 1, the default genetic code
    in scikit-bio):

    >>> from skbio import GeneticCode
    >>> GeneticCode.from_ncbi()
    GeneticCode (Standard)
    -------------------------------------------------------------------------
      AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Starts = ---M---------------M---------------M----------------------------
    Base1  = UUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2  = UUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGG
    Base3  = UCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAG

    Get a different NCBI genetic code (25):

    >>> GeneticCode.from_ncbi(25)
    GeneticCode (Candidate Division SR1 and Gracilibacteria)
    -------------------------------------------------------------------------
      AAs  = FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Starts = ---M-------------------------------M---------------M------------
    Base1  = UUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2  = UUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGG
    Base3  = UCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAG

    Define a custom genetic code:

    >>> GeneticCode('M' * 64, '-' * 64)
    GeneticCode
    -------------------------------------------------------------------------
      AAs  = MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    Starts = ----------------------------------------------------------------
    Base1  = UUUUUUUUUUUUUUUUCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2  = UUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGGUUUUCCCCAAAAGGGG
    Base3  = UCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAGUCAG

    Translate an RNA sequence to protein using NCBI's standard genetic code:

    >>> from skbio import RNA
    >>> rna = RNA('AUGCCACUUUAA')
    >>> GeneticCode.from_ncbi().translate(rna)
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

    """

    _num_codons = 64
    _radix_multiplier = np.asarray([16, 4, 1], dtype=np.uint8)
    _start_stop_options = ["ignore", "optional", "require"]
    __offset_table = None

    @classproperty
    def _offset_table(cls):
        if cls.__offset_table is None:
            # create lookup table that is filled with 255 everywhere except for
            # indices corresponding to U, C, A, and G. 255 was chosen to
            # represent invalid character offsets because it will create an
            # invalid (out of bounds) index into `amino_acids` which should
            # error noisily. this is important in case the valid definite
            # IUPAC RNA characters change in the future and the assumptions
            # currently made by the code become invalid
            table = np.empty(ord(b"U") + 1, dtype=np.uint8)
            table.fill(255)
            table[ord(b"U")] = 0
            table[ord(b"C")] = 1
            table[ord(b"A")] = 2
            table[ord(b"G")] = 3
            cls.__offset_table = table
        return cls.__offset_table

    @classonlymethod
    def from_ncbi(cls, table_id=1):
        """Return NCBI genetic code specified by table ID.

        Parameters
        ----------
        table_id : int, optional
            Table ID of the NCBI genetic code to return.

        Returns
        -------
        GeneticCode
            NCBI genetic code specified by `table_id`.

        Notes
        -----
        The table IDs and genetic codes available in this method and used
        throughout the examples are defined in [1]_.

        References
        ----------
        .. [1] http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

        Examples
        --------
        Get the NCBI thraustochytrium mitochondrial genetic code (23):

        >>> tmgc = GeneticCode.from_ncbi(23)
        >>> tmgc.name
        'Thraustochytrium Mitochondrial'

        """
        if table_id not in _ncbi_genetic_codes:
            raise ValueError(
                "`table_id` must be one of %r, not %r"
                % (sorted(_ncbi_genetic_codes), table_id)
            )
        return _ncbi_genetic_codes[table_id]

    @classproperty
    def reading_frames(cls):
        """Six possible reading frames.

        Reading frames are ordered:

        * 1 (forward)
        * 2 (forward)
        * 3 (forward)
        * -1 (reverse)
        * -2 (reverse)
        * -3 (reverse)

        This property can be passed into
        ``GeneticCode.translate(reading_frame)``.

        Returns
        -------
        list (int)
            Six possible reading frames.

        """
        return [1, 2, 3, -1, -2, -3]

    @property
    def name(self):
        """Genetic code name.

        This is simply metadata and does not affect the functionality of the
        genetic code itself.

        Returns
        -------
        str
            Genetic code name.

        """
        return self._name

    def __init__(self, amino_acids, starts, name=""):
        self._set_amino_acids(amino_acids)
        self._set_starts(starts)
        self._name = name

    def _set_amino_acids(self, amino_acids):
        amino_acids = Protein(amino_acids)

        if len(amino_acids) != self._num_codons:
            raise ValueError(
                "`amino_acids` must be length %d, not %d"
                % (self._num_codons, len(amino_acids))
            )
        indices = (amino_acids.values == b"M").nonzero()[0]
        if indices.size < 1:
            raise ValueError(
                "`amino_acids` must contain at least one M " "(methionine) character"
            )
        self._amino_acids = amino_acids
        self._m_character_codon = self._index_to_codon(indices[0])

    def _set_starts(self, starts):
        starts = Protein(starts)

        if len(starts) != self._num_codons:
            raise ValueError(
                "`starts` must be length %d, not %d" % (self._num_codons, len(starts))
            )
        if (starts.values == b"M").sum() + (starts.values == b"-").sum() != len(starts):
            # to prevent the user from accidentally swapping `starts` and
            # `amino_acids` and getting a translation back
            raise ValueError("`starts` may only contain M and - characters")

        self._starts = starts

        indices = (self._starts.values == b"M").nonzero()[0]
        codons = np.empty((indices.size, 3), dtype=np.uint8)
        for i, index in enumerate(indices):
            codons[i] = self._index_to_codon(index)
        self._start_codons = codons

    def _index_to_codon(self, index):
        """Convert AA index (0-63) to codon encoded in offsets (0-3)."""
        codon = np.empty(3, dtype=np.uint8)
        for i, multiplier in enumerate(self._radix_multiplier):
            offset, index = divmod(index, multiplier)
            codon[i] = offset
        return codon

    def __str__(self):
        """Return string representation of the genetic code.

        Returns
        -------
        str
            Genetic code in NCBI genetic code format.

        Notes
        -----
        Representation uses NCBI genetic code format defined in [1]_.

        References
        ----------
        .. [1] http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

        """
        return self._build_repr(include_name=False)

    def __repr__(self):
        """Return string representation of the genetic code.

        Returns
        -------
        str
            Genetic code in NCBI genetic code format.

        Notes
        -----
        Representation uses NCBI genetic code format defined in [1]_ preceded
        by a header. If the genetic code has a name, it will be included in the
        header.

        References
        ----------
        .. [1] http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

        """
        return self._build_repr(include_name=True)

    def _build_repr(self, include_name):
        lines = ElasticLines()

        if include_name:
            name_line = self.__class__.__name__
            if len(self.name) > 0:
                name_line += " (%s)" % self.name
            lines.add_line(name_line)
            lines.add_separator()

        lines.add_line("  AAs  = %s" % str(self._amino_acids))
        lines.add_line("Starts = %s" % str(self._starts))
        base1 = "U" * 16 + "C" * 16 + "A" * 16 + "G" * 16
        lines.add_line("Base1  = %s" % base1)
        base2 = ("U" * 4 + "C" * 4 + "A" * 4 + "G" * 4) * 4
        lines.add_line("Base2  = %s" % base2)
        base3 = "UCAG" * 16
        lines.add_line("Base3  = %s" % base3)

        return lines.to_str()

    def __eq__(self, other):
        """Determine if the genetic code is equal to another.

        Genetic codes are equal if they are *exactly* the same type and
        defined by the same `amino_acids` and `starts`. A genetic code's name
        (accessed via ``name`` property) does not affect equality.

        Parameters
        ----------
        other : GeneticCode
            Genetic code to test for equality against.

        Returns
        -------
        bool
            Indicates whether the genetic code is equal to `other`.

        Examples
        --------
        NCBI genetic codes 1 and 2 are not equal:

        >>> GeneticCode.from_ncbi(1) == GeneticCode.from_ncbi(2)
        False

        Define a custom genetic code:

        >>> gc = GeneticCode('M' * 64, '-' * 64)

        Define a second genetic code with the same `amino_acids` and `starts`.
        Note that the presence of a name does not make the genetic codes
        unequal:

        >>> named_gc = GeneticCode('M' * 64, '-' * 64, name='example name')
        >>> gc == named_gc
        True

        """
        if self.__class__ != other.__class__:
            return False
        # convert Protein to str so that metadata is ignored in comparison. we
        # only care about the sequence data defining the genetic code
        if str(self._amino_acids) != str(other._amino_acids):
            return False
        if str(self._starts) != str(other._starts):
            return False
        return True

    def __ne__(self, other):
        """Determine if the genetic code is not equal to another.

        Genetic codes are not equal if their type, `amino_acids`, or `starts`
        differ. A genetic code's name (accessed via ``name`` property) does not
        affect equality.

        Parameters
        ----------
        other : GeneticCode
            Genetic code to test for inequality against.

        Returns
        -------
        bool
            Indicates whether the genetic code is not equal to `other`.

        """
        return not (self == other)

    def translate(self, sequence, reading_frame=1, start="ignore", stop="ignore"):
        """Translate RNA sequence into protein sequence.

        Parameters
        ----------
        sequence : RNA
            RNA sequence to translate.
        reading_frame : {1, 2, 3, -1, -2, -3}
            Reading frame to use in translation. 1, 2, and 3 are forward frames
            and -1, -2, and -3 are reverse frames. If reverse (negative), will
            reverse complement the sequence before translation.
        start : {'ignore', 'require', 'optional'}
            How to handle start codons:

            * "ignore": translation will start from the beginning of the
              reading frame, regardless of the presence of a start codon.

            * "require": translation will start at the first start codon in
              the reading frame, ignoring all prior positions. The first amino
              acid in the translated sequence will *always* be methionine
              (M character), even if an alternative start codon was used in
              translation. This behavior most closely matches the underlying
              biology since fMet doesn't have a corresponding IUPAC character.
              If a start codon does not exist, a ``ValueError`` is raised.

            * "optional": if a start codon exists in the reading frame, matches
              the behavior of "require". If a start codon does not exist,
              matches the behavior of "ignore".

        stop : {'ignore', 'require', 'optional'}
            How to handle stop codons:

            * "ignore": translation will ignore the presence of stop codons and
              translate to the end of the reading frame.

            * "require": translation will terminate at the first stop codon.
              The stop codon will not be included in the translated sequence.
              If a stop codon does not exist, a ``ValueError`` is raised.

            * "optional": if a stop codon exists in the reading frame, matches
              the behavior of "require". If a stop codon does not exist,
              matches the behavior of "ignore".

        Returns
        -------
        Protein
            Translated sequence.

        See Also
        --------
        translate_six_frames

        Notes
        -----
        Input RNA sequence metadata are included in the translated protein
        sequence. Positional metadata are not included.

        Examples
        --------
        Translate RNA into protein using NCBI's standard genetic code (table ID
        1, the default genetic code in scikit-bio):

        >>> from skbio import RNA, GeneticCode
        >>> rna = RNA('AGUAUUCUGCCACUGUAAGAA')
        >>> sgc = GeneticCode.from_ncbi()
        >>> sgc.translate(rna)
        Protein
        --------------------------
        Stats:
            length: 7
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: True
        --------------------------
        0 SILPL*E

        In this command, we used the default ``start`` behavior, which starts
        translation at the beginning of the reading frame, regardless of the
        presence of a start codon. If we specify "require", translation will
        start at the first start codon in the reading frame (in this example,
        CUG), ignoring all prior positions:

        >>> sgc.translate(rna, start='require')
        Protein
        --------------------------
        Stats:
            length: 5
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: True
        --------------------------
        0 MPL*E

        Note that the codon coding for L (CUG) is an alternative start codon in
        this genetic code. Since we specified "require" mode, methionine (M)
        was used in place of the alternative start codon (L). This behavior
        most closely matches the underlying biology since fMet doesn't have a
        corresponding IUPAC character.

        Translate the same RNA sequence, also specifying that translation
        terminate at the first stop codon in the reading frame:

        >>> sgc.translate(rna, start='require', stop='require')
        Protein
        --------------------------
        Stats:
            length: 3
            has gaps: False
            has degenerates: False
            has definites: True
            has stops: False
        --------------------------
        0 MPL

        Passing "require" to both ``start`` and ``stop`` trims the translation
        to the CDS (and in fact requires that one is present in the reading
        frame). Changing the reading frame to 2 causes an exception to be
        raised because a start codon doesn't exist in the reading frame:

        >>> sgc.translate(rna, start='require', stop='require',
        ...               reading_frame=2) # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ValueError: ...

        """
        self._validate_translate_inputs(sequence, reading_frame, start, stop)

        offset = abs(reading_frame) - 1
        if reading_frame < 0:
            sequence = sequence.reverse_complement()

        # Translation strategy:
        #
        #   1. Obtain view of underlying sequence bytes from the beginning of
        #      the reading frame.
        #   2. Convert bytes to offsets (0-3, base 4 since there are only 4
        #      characters allowed: UCAG).
        #   3. Reshape byte vector into (N, 3), where N is the number of codons
        #      in the reading frame. Each row represents a codon in the
        #      sequence.
        #   4. (Optional) Find start codon in the reading frame and trim to
        #      this position. Replace start codon with M codon.
        #   5. Convert each codon (encoded as offsets) into an index
        #      corresponding to an amino acid (0-63).
        #   6. Obtain translated sequence by indexing into the amino acids
        #      vector (`amino_acids`) using the indices defined in step 5.
        #   7. (Optional) Find first stop codon and trim to this position.
        data = sequence.values[offset:].view(np.uint8)
        # since advanced indexing is used with an integer ndarray, a copy is
        # always returned. thus, the in-place modification made below
        # (replacing the start codon) is safe.
        data = self._offset_table[data]
        data = data[: data.size // 3 * 3].reshape((-1, 3))

        if start in {"require", "optional"}:
            start_codon_index = data.shape[0]
            for start_codon in self._start_codons:
                indices = np.all(data == start_codon, axis=1).nonzero()[0]

                if indices.size > 0:
                    first_index = indices[0]
                    if first_index < start_codon_index:
                        start_codon_index = first_index

            if start_codon_index != data.shape[0]:
                data = data[start_codon_index:]
                data[0] = self._m_character_codon
            elif start == "require":
                self._raise_require_error("start", reading_frame)

        indices = (data * self._radix_multiplier).sum(axis=1)
        translated = self._amino_acids.values[indices]

        if stop in {"require", "optional"}:
            stop_codon_indices = (translated == b"*").nonzero()[0]
            if stop_codon_indices.size > 0:
                translated = translated[: stop_codon_indices[0]]
            elif stop == "require":
                self._raise_require_error("stop", reading_frame)

        metadata = None
        if sequence.has_metadata():
            metadata = sequence.metadata

        # turn off validation because `translated` is guaranteed to be valid
        return Protein(translated, metadata=metadata, validate=False)

    def _validate_translate_inputs(self, sequence, reading_frame, start, stop):
        if not isinstance(sequence, RNA):
            raise TypeError(
                "Sequence to translate must be RNA, not %s" % type(sequence).__name__
            )

        if reading_frame not in self.reading_frames:
            raise ValueError(
                "`reading_frame` must be one of %r, not %r"
                % (self.reading_frames, reading_frame)
            )

        for name, value in ("start", start), ("stop", stop):
            if value not in self._start_stop_options:
                raise ValueError(
                    "`%s` must be one of %r, not %r"
                    % (name, self._start_stop_options, value)
                )

        if sequence.has_gaps():
            raise ValueError(
                "scikit-bio does not support translation of " "gapped sequences."
            )

        if sequence.has_degenerates():
            raise NotImplementedError(
                "scikit-bio does not currently support "
                "translation of degenerate sequences."
                "`RNA.expand_degenerates` can be used "
                "to obtain all definite versions "
                "of a degenerate sequence."
            )

    def _raise_require_error(self, name, reading_frame):
        raise ValueError(
            "Sequence does not contain a %s codon in the "
            "current reading frame (`reading_frame=%d`). Presence "
            "of a %s codon is required with `%s='require'`"
            % (name, reading_frame, name, name)
        )

    def translate_six_frames(self, sequence, start="ignore", stop="ignore"):
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
        sequence : RNA
            RNA sequence to translate.
        start : {'ignore', 'require', 'optional'}
            How to handle start codons. See ``GeneticCode.translate`` for
            details.
        stop : {'ignore', 'require', 'optional'}
            How to handle stop codons. See ``GeneticCode.translate`` for
            details.

        Yields
        ------
        Protein
            Translated sequence in the current reading frame.

        See Also
        --------
        translate

        Notes
        -----
        This method is faster than (and equivalent to) performing six
        independent translations using, for example:

        ``(gc.translate(seq, reading_frame=rf)
        for rf in GeneticCode.reading_frames)``

        Input RNA sequence metadata are included in each translated protein
        sequence. Positional metadata are not included.

        Examples
        --------
        Translate RNA into protein using the six possible reading frames and
        NCBI's standard genetic code (table ID 1, the default genetic code in
        scikit-bio):

        >>> from skbio import RNA, GeneticCode
        >>> rna = RNA('AUGCCACUUUAA')
        >>> sgc = GeneticCode.from_ncbi()
        >>> for protein in sgc.translate_six_frames(rna):
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
        rc = sequence.reverse_complement()

        for reading_frame in range(1, 4):
            yield self.translate(
                sequence, reading_frame=reading_frame, start=start, stop=stop
            )
        for reading_frame in range(1, 4):
            yield self.translate(
                rc, reading_frame=reading_frame, start=start, stop=stop
            )


# defined at http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
_ncbi_genetic_codes = {
    1: GeneticCode(
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M---------------M---------------M----------------------------",
        "Standard",
    ),
    2: GeneticCode(
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "--------------------------------MMMM---------------M------------",
        "Vertebrate Mitochondrial",
    ),
    3: GeneticCode(
        "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "----------------------------------MM----------------------------",
        "Yeast Mitochondrial",
    ),
    4: GeneticCode(
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM---------------M------------MMMM---------------M------------",
        "Mold, Protozoan, and Coelenterate Mitochondrial, and "
        "Mycoplasma/Spiroplasma",
    ),
    5: GeneticCode(
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "---M----------------------------MMMM---------------M------------",
        "Invertebrate Mitochondrial",
    ),
    6: GeneticCode(
        "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "Ciliate, Dasycladacean and Hexamita Nuclear",
    ),
    9: GeneticCode(
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "-----------------------------------M---------------M------------",
        "Echinoderm and Flatworm Mitochondrial",
    ),
    10: GeneticCode(
        "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "Euplotid Nuclear",
    ),
    11: GeneticCode(
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M---------------M------------MMMM---------------M------------",
        "Bacterial, Archaeal and Plant Plastid",
    ),
    12: GeneticCode(
        "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-------------------M---------------M----------------------------",
        "Alternative Yeast Nuclear",
    ),
    13: GeneticCode(
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        "---M------------------------------MM---------------M------------",
        "Ascidian Mitochondrial",
    ),
    14: GeneticCode(
        "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "Alternative Flatworm Mitochondrial",
    ),
    16: GeneticCode(
        "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "Chlorophycean Mitochondrial",
    ),
    21: GeneticCode(
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "-----------------------------------M---------------M------------",
        "Trematode Mitochondrial",
    ),
    22: GeneticCode(
        "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "Scenedesmus obliquus Mitochondrial",
    ),
    23: GeneticCode(
        "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--------------------------------M--M---------------M------------",
        "Thraustochytrium Mitochondrial",
    ),
    24: GeneticCode(
        "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        "---M---------------M---------------M---------------M------------",
        "Pterobranchia Mitochondrial",
    ),
    25: GeneticCode(
        "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M-------------------------------M---------------M------------",
        "Candidate Division SR1 and Gracilibacteria",
    ),
}
