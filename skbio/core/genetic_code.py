#!/usr/bin/env python
r"""
Genetic Code (:mod:`skbio.core.genetic_code`)
=============================================

.. currentmodule:: skbio.core.genetic_code

This module defines the ``GeneticCode`` class, which represents an immutable
object that translates RNA or DNA strings to amino acid sequences.


Classes
-------

.. autosummary::
   :toctree: generated/

   GeneticCode

Examples
--------

Creating and using a ``GeneticCode`` object

>>> from skbio.core.genetic_code import GeneticCodes
>>> from pprint import pprint
>>> sgc = GeneticCodes[1]
>>> sgc
GeneticCode(FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG)
>>> sgc['UUU'] == 'F'
True
>>> sgc['TTT'] == 'F'
True
>>> sgc['F'] == ['TTT', 'TTC']          #in arbitrary order
True
>>> sgc['*'] == ['TAA', 'TAG', 'TGA']   #in arbitrary order
True

Retrieving the anticodons of the object

>>> pprint(sgc.anticodons)
{'*': ['TTA', 'CTA', 'TCA'],
 'A': ['AGC', 'GGC', 'TGC', 'CGC'],
 'C': ['ACA', 'GCA'],
 'D': ['ATC', 'GTC'],
 'E': ['TTC', 'CTC'],
 'F': ['AAA', 'GAA'],
 'G': ['ACC', 'GCC', 'TCC', 'CCC'],
 'H': ['ATG', 'GTG'],
 'I': ['AAT', 'GAT', 'TAT'],
 'K': ['TTT', 'CTT'],
 'L': ['TAA', 'CAA', 'AAG', 'GAG', 'TAG', 'CAG'],
 'M': ['CAT'],
 'N': ['ATT', 'GTT'],
 'P': ['AGG', 'GGG', 'TGG', 'CGG'],
 'Q': ['TTG', 'CTG'],
 'R': ['ACG', 'GCG', 'TCG', 'CCG', 'TCT', 'CCT'],
 'S': ['AGA', 'GGA', 'TGA', 'CGA', 'ACT', 'GCT'],
 'T': ['AGT', 'GGT', 'TGT', 'CGT'],
 'V': ['AAC', 'GAC', 'TAC', 'CAC'],
 'W': ['CCA'],
 'Y': ['ATA', 'GTA']}

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
import re

from collections import defaultdict

from skbio.core.exception import GeneticCodeInitError, InvalidCodonError
from skbio.core.sequence import Protein

# py3k compatibility
try:
    from string import maketrans
except ImportError:
    maketrans = str.maketrans

_dna_trans = maketrans('TCAG', 'AGTC')


def _simple_rc(seq):
    """simple reverse-complement: works only on unambiguous uppercase DNA"""
    return seq.translate(_dna_trans)[::-1]


class GeneticCode(object):

    """Class to hold codon to amino acid mapping, and vice versa.

    Attributes
    ----------
    code_sequence
    id
    name
    start_codon_sequence
    start_codons
    codons
    synonyms
    sense_codons
    anticodons
    blocks

    Parameters
    ----------
    code_sequence : str
        64-character string containing NCBI representation.
    id : str, optional
        identifier for the object.
    name : str, optional
        name for the object.
    start_codon_sequence : str, optional
        starting point for the codon sequence.

    Returns
    -------
    GeneticCode
        initialized ``GeneticCode`` object.

    Raises
    ------
    GeneticCodeInitError
        If the length of `code_sequence` is different to `64`.

    Methods
    -------
    changes
    get_stop_indices
    is_start
    is_stop
    translate_six_frames
    translate
    __repr__
    __getitem__
    __str__
    __eq__

    Examples
    --------
    >>> from skbio.core.genetic_code import GeneticCode
    >>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSR'
    ...                   'RVVVVAAAADDEEGGGG')

    .. note:: `*` is used to denote termination as per the NCBI standard.
        Although the genetic code objecs convert DNA to RNA and vice versa,
        lists of codons that they produce will be provided in DNA format.

    """
    # class data: need the bases, the list of codons in UUU -> GGG order, and
    # a mapping from positions in the list back to codons. These should be the
    # same for all GeneticCode instances, and are immutable (therefore
    # private).
    _codons = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]

    def __init__(self, code_sequence, id=None, name=None,
                 start_codon_sequence=None):
        if len(code_sequence) != 64:
            raise GeneticCodeInitError("code_sequence: %s has length %d, but "
                                       "expected 64" % (code_sequence,
                                                        len(code_sequence)))

        self.code_sequence = code_sequence
        self.id = id
        self.name = name
        self.start_codon_sequence = start_codon_sequence
        start_codons = {}
        if start_codon_sequence is not None:
            for codon, aa in zip(self._codons, start_codon_sequence):
                if aa != '-':
                    start_codons[codon] = aa
        self.start_codons = start_codons
        codon_lookup = {key: value for (key, value) in zip(self._codons,
                                                           code_sequence)}
        self.codons = codon_lookup

        # create synonyms for each aa
        aa_lookup = defaultdict(list)
        for codon in self._codons:
            aa = codon_lookup[codon]
            aa_lookup[aa].append(codon)
        self.synonyms = dict(aa_lookup)
        sense_codons = codon_lookup.copy()

        # create sense codons
        stop_codons = self['*']
        for c in stop_codons:
            del sense_codons[c]
        self.sense_codons = sense_codons

        # create anticodons
        ac = {}
        for aa, codons in self.synonyms.items():
            ac[aa] = [_simple_rc(element) for element in codons]
        self.anticodons = ac

    def _analyze_quartet(self, codons, aa):
        """Analyzes a quartet of codons and amino acids: returns list of lists.

        Each list contains one block, splitting at purine/pyrimidine boundary
        if necessary.

        codons should be a list of 4 codons.
        aa should be a list of 4 amino acid symbols.

        Possible states:
            - All amino acids are the same: returns list of one quartet.
            - Two groups of 2 aa: returns list of two doublets.
            - One group of 2 and 2 groups of 1: list of one doublet, 2 singles.
            - 4 groups of 1: four singles.

        Note: codon blocks like Ile in the standard code (AUU, AUC, AUA) will
        be split when they cross the R/Y boundary, so [[AUU, AUC], [AUA]]. This
        would also apply to a block like AUC AUA AUG -> [[AUC],[AUA,AUG]],
        although this latter pattern is not observed in the standard code.
        """
        if aa[0] == aa[1]:
            first_doublet = True
        else:
            first_doublet = False
        if aa[2] == aa[3]:
            second_doublet = True
        else:
            second_doublet = False
        if first_doublet and second_doublet and aa[1] == aa[2]:
            return [codons]
        else:
            blocks = []
            if first_doublet:
                blocks.append(codons[:2])
            else:
                blocks.extend([[codons[0]], [codons[1]]])
            if second_doublet:
                blocks.append(codons[2:])
            else:
                blocks.extend([[codons[2]], [codons[3]]])
            return blocks

    def _get_blocks(self):
        """Returns list of lists of codon blocks in the genetic code.

        A codon block can be:
            - a quartet, if all 4 XYn codons have the same amino acid.
            - a doublet, if XYt and XYc or XYa and XYg have the same aa.
            - a singlet, otherwise.

        Returns
        -------
        list
            Returns a list of the quartets, doublets, and singlets in the order
            UUU -> GGG.

        Notes
        -----
        A doublet cannot span the purine/pyrimidine boundary, and a quartet
        cannot span the boundary between two codon blocks whose first two bases
        differ.

        """
        if hasattr(self, '_blocks'):
            return self._blocks
        else:
            blocks = []
            curr_codons = []
            curr_aa = []
            for index, codon, aa in zip(range(64), self._codons,
                                        self.code_sequence):
                # we're in a new block if it's a new quartet or a different aa
                new_quartet = not index % 4
                if new_quartet and curr_codons:
                    blocks.extend(self._analyze_quartet(curr_codons, curr_aa))
                    curr_codons = []
                    curr_aa = []
                curr_codons.append(codon)
                curr_aa.append(aa)
            # don't forget to append last block
            if curr_codons:
                blocks.extend(self._analyze_quartet(curr_codons, curr_aa))
            self._blocks = blocks
            return self._blocks

    blocks = property(_get_blocks)

    def __str__(self):
        """Returns code_sequence that constructs the GeneticCode

        .. shownumpydoc
        """
        return self.code_sequence

    def __repr__(self):
        """Returns reconstructable representation of the GeneticCode

        .. shownumpydoc
        """
        return 'GeneticCode(%s)' % str(self)

    def __eq__(self, other):
        """ Allows two GeneticCode objects to be compared to each other.

        Two GeneticCode objects are equal if they have equal code_sequences.

        .. shownumpydoc
        """
        if not isinstance(other, GeneticCode):
            return False
        return self.code_sequence == other.code_sequence

    def __ne__(self, other):
        """Required in Py2."""
        return not self == other

    def __getitem__(self, item):
        """Returns amino acid corresponding to codon, or codons for an aa.

        Returns [] for empty list of codons, 'X' for unknown amino acid.

        .. shownumpydoc
        """
        item = str(item)
        if len(item) == 1:  # amino acid
            return self.synonyms.get(item, [])
        elif len(item) == 3:  # codon
            key = item.upper()
            key = key.replace('U', 'T')
            return self.codons.get(key, 'X')
        else:
            raise InvalidCodonError("Codon or aa %s has wrong length" % item)

    def translate(self, nucleotide_sequence, start=0):
        """Translate nucleotide to protein sequence

        Parameters
        ----------
        nucleotide_sequence : NucleotideSequence
            sequence to be translated
        start : int, optional
            position to begin translation

        Returns
        -------
        ProteinSequence
            translation of nucleotide_sequence

        Notes
        -----
        ``translate`` returns the translation of the entire sequence, (i.e., of
        ``nucleotide_sequence[start:]``). It is the user's responsibility to
        trim to an open reading frame, either from the input or using the
        output, if that is desired.

        See Also
        --------
        translate_six_frames

        Examples
        --------
        >>> from skbio.core.genetic_code import GeneticCode
        >>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS'
        ...                   'RRVVVVAAAADDEEGGGG')
        >>> sgc.translate('AUGCAUGACUUUUGA', 1)
        <ProteinSequence: CMTF (length: 4)>

        """
        if len(nucleotide_sequence) == 0:
            return Protein('')
        if start + 1 > len(nucleotide_sequence):
            raise ValueError("Translation starts after end of"
                             "NucleotideSequence")

        translation = []
        for i in range(start, len(nucleotide_sequence) - 2, 3):
            translation.append(self[nucleotide_sequence[i:i + 3]])
        translation = Protein(''.join(translation))

        return translation

    def get_stop_indices(self, nucleotide_sequence, start=0):
        """returns indexes for stop codons in the specified frame

        Parameters
        ----------
        nucleotide_sequence : str, NucleotideSequence
            sequence to be scanned for stop codons
        start : int, optional
            position where the search begins.

        Returns
        -------
        list
            indices of the stop codons.

        Examples
        --------
        >>> from skbio.core.genetic_code import GeneticCode
        >>> from skbio.core.sequence import DNA
        >>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS'
        ...                   'RRVVVVAAAADDEEGGGG')
        >>> seq = DNA('ATGCTAACATAAA')
        >>> sgc.get_stop_indices(seq, 0)
        [9]

        """
        stops = self['*']
        stop_pattern = '(%s)' % '|'.join(stops)
        stop_pattern = re.compile(stop_pattern)
        seq = str(nucleotide_sequence)
        found = [hit.start() for hit in stop_pattern.finditer(seq)]
        found = [index for index in found if index % 3 == start]
        return found

    def translate_six_frames(self, nucleotide_sequence):
        """Translate nucleotide to protein sequences for all six reading frames

        Parameters
        ----------
        nucleotide_sequence : NucleotideSequence
            sequence to be scanned for stop codons

        Returns
        -------
        list
            the six translated ProteinSequence objects

        See Also
        --------
        translate

        Examples
        --------
        >>> from skbio.core.genetic_code import GeneticCode
        >>> from skbio.core.sequence import RNA
        >>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS'
        ...                   'RRVVVVAAAADDEEGGGG')
        >>> results = sgc.translate_six_frames(RNA('AUGCUAACAUAAA'))
        >>> for e in results: e
        <ProteinSequence: MLT* (length: 4)>
        <ProteinSequence: C*HK (length: 4)>
        <ProteinSequence: ANI (length: 3)>
        <ProteinSequence: FMLA (length: 4)>
        <ProteinSequence: LC*H (length: 4)>
        <ProteinSequence: YVS (length: 3)>

        """
        rc_nucleotide_sequence = nucleotide_sequence.rc()
        results = []
        for start in range(3):
            translation = self.translate(nucleotide_sequence, start)
            results.append(translation)

        for start in range(3):
            translation = self.translate(rc_nucleotide_sequence, start)
            results.append(translation)

        return results

    def is_start(self, codon):
        """Checks if codon is a start codon

        Parameters
        ----------
        codon : str
            codon string

        Returns
        -------
        bool
            ``True`` if codon is a start codon, ``False`` otherwise

        Examples
        --------
        >>> from skbio.core.genetic_code import GeneticCode
        >>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS'
        ...                   'RRVVVVAAAADDEEGGGG')
        >>> sgc.is_start('ATG')
        False
        >>> sgc.is_start('AAA')
        False

        """
        fixed_codon = codon.upper().replace('U', 'T')
        return fixed_codon in self.start_codons

    def is_stop(self, codon):
        """Checks if codon is a stop codon

        Parameters
        ----------
        codon : str
            codon string

        Returns
        -------
        bool
            ``True`` if codon is a stop codon, ``False`` otherwise

        Examples
        --------
        >>> from skbio.core.genetic_code import GeneticCode
        >>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS'
        ...                   'RRVVVVAAAADDEEGGGG')
        >>> sgc.is_stop('UAA')
        True
        >>> sgc.is_stop('AAA')
        False

        """
        return self[codon] == '*'

    def changes(self, other):
        """Returns dictionary of coddons that differ

        Parameters
        ----------
        other : GeneticCode
           genetic code object

        Returns
        -------
        dict
            Returns a dictionary of the form ``{codon:'XY'}`` for codons that
            differ. X is the string representation of the amino acid in the
            object calling this method, Y is the string representation of the
            amino acid in `other`. Always returns a 2-character string.

        Examples
        --------
        >>> from skbio.core.genetic_code import GeneticCode
        >>> from pprint import pprint
        >>> sgc = GeneticCode('FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS*'
        ...                   '*VVVVAAAADDEEGGGG')
        >>> pprint(sgc.changes('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTT'
        ...                    'TNNKKSSRRVVVVAAAADDEEGGGG'))
        {'AGA': '*R', 'AGG': '*R', 'ATA': 'MI', 'TGA': 'W*'}

        """
        changes = {}
        try:
            other_code = other.code_sequence
        except AttributeError:  # try using other directly as sequence
            other_code = other
        for codon, old, new in zip(self._codons, self.code_sequence,
                                   other_code):
            if old != new:
                changes[codon] = old + new
        return changes


NCBIGeneticCodeData = [GeneticCode(*data) for data in [
    [
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        1,
        'Standard Nuclear',
        '---M---------------M---------------M----------------------------',
    ],
    [
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
        2,
        'Vertebrate Mitochondrial',
        '--------------------------------MMMM---------------M------------',
    ],
    [
        'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        3,
        'Yeast Mitochondrial',
        '----------------------------------MM----------------------------',
    ],
    [
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        4,
        'Mold, Protozoan, and Coelenterate Mitochondrial, and Mycoplasma/'
        'Spiroplasma Nuclear',
        '--MM---------------M------------MMMM---------------M------------',
    ],
    [
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
        5,
        'Invertebrate Mitochondrial',
        '---M----------------------------MMMM---------------M------------',
    ],
    [
        'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        6,
        'Ciliate, Dasycladacean and Hexamita Nuclear',
        '-----------------------------------M----------------------------',
    ],
    [
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        9,
        'Echinoderm and Flatworm Mitochondrial',
        '-----------------------------------M---------------M------------',
    ],
    [
        'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        10,
        'Euplotid Nuclear',
        '-----------------------------------M----------------------------',
    ],
    [
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        11,
        'Bacterial Nuclear and Plant Plastid',
        '---M---------------M------------MMMM---------------M------------',
    ],
    [
        'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        12,
        'Alternative Yeast Nuclear',
        '-------------------M---------------M----------------------------',
    ],
    [
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
        13,
        'Ascidian Mitochondrial',
        '-----------------------------------M----------------------------',
    ],
    [
        'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        14,
        'Alternative Flatworm Mitochondrial',
        '-----------------------------------M----------------------------',
    ],
    [
        'FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        15,
        'Blepharisma Nuclear',
        '-----------------------------------M----------------------------',
    ],
    [
        'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        16,
        'Chlorophycean Mitochondrial',
        '-----------------------------------M----------------------------',
    ],
    [
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        20,
        'Trematode Mitochondrial',
        '-----------------------------------M---------------M------------',
    ],
    [
        'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        22,
        'Scenedesmus obliquus Mitochondrial',
        '-----------------------------------M----------------------------',
    ],
    [
        'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        23,
        'Thraustochytrium Mitochondrial',
    ],
]]

# build dict of GeneticCodes keyed by ID (as int, not str)
GeneticCodes = dict([(i.id, i) for i in NCBIGeneticCodeData])
# add str versions for convenience
for key, value in list(GeneticCodes.items()):
    GeneticCodes[str(key)] = value

DEFAULT = GeneticCodes[1]
