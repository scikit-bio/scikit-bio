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

>>> from skbio.core.genetic_code import GeneticCode
>>> from pprint import pprint
>>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAA'
...                   'ADDEEGGGG')
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

from skbio.core.exception import GeneticCodeInitError, InvalidCodonError

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
    sixframes
    translate
    __repr__
    __getitem__
    __cmp__
    __str__


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
        if start_codon_sequence:
            for codon, aa in zip(self._codons, start_codon_sequence):
                if aa != '-':
                    start_codons[codon] = aa
        self.start_codons = start_codons
        codon_lookup = dict(zip(self._codons, code_sequence))
        self.codons = codon_lookup
        # create synonyms for each aa
        aa_lookup = {}
        for codon in self._codons:
            aa = codon_lookup[codon]
            if aa not in aa_lookup:
                aa_lookup[aa] = [codon]
            else:
                aa_lookup[aa].append(codon)
        self.synonyms = aa_lookup
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

        Each list contains one block, splitting at R/Y if necessary.

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

    def __cmp__(self, other):
        """ Allows two GeneticCode objects to be compared to each other.

        Two GeneticCode objects are equal if they have equal code_sequences.

        .. shownumpydoc
        """
        return cmp(str(self), str(other))

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

    def translate(self, dna, start=0):
        """Translates DNA to protein with current GeneticCode.

        Translates the entire sequence: it is the caller's responsibility to
        find open reading frames.

        Parameters
        ----------
        dna : str
            a string of nucleotides
        start : int, optional
            position to begin translation (used to implement frames)

        Returns
        -------
        str
            string containing amino acid sequence.

        Examples
        --------
        >>> from skbio.core.genetic_code import GeneticCode
        >>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS'
        ...                   'RRVVVVAAAADDEEGGGG')
        >>> sgc.translate('AUGCAUGACUUUUGA', 1)
        'CMTF'

        """
        # NOTE: should return Protein object when we have a class for it.
        if not dna:
            return ''
        if start + 1 > len(dna):
            raise ValueError("Translation starts after end of RNA")
        return ''.join([self[dna[i:i + 3]] for i in
                        range(start, len(dna) - 2, 3)])

    def get_stop_indices(self, dna, start=0):
        """returns indexes for stop codons in the specified frame

        Parameters
        ----------
        dna : str
            a string of nucleotides.
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
        seq = str(dna)
        found = [hit.start() for hit in stop_pattern.finditer(seq)]
        found = [index for index in found if index % 3 == start]
        return found

    def sixframes(self, dna):
        """Returns six-frame translation as a dictionary object

        Parameters
        ----------
        dna : str
            a string of nucleotides.

        Returns
        -------
        dict
            dictionary of where the keys are the frames and the values are the
            translations i. e. ``{frame:translation}``.

        Examples
        --------
        >>> from skbio.core.genetic_code import GeneticCode
        >>> from skbio.core.sequence import RNA
        >>> sgc = GeneticCode('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS'
        ...                   'RRVVVVAAAADDEEGGGG')
        >>> sgc.sixframes(RNA('AUGCUAACAUAAA'))
        ['MLT*', 'C*HK', 'ANI', 'FMLA', 'LC*H', 'YVS']

        """
        reverse = dna.rc()
        return [self.translate(dna, start) for start in range(3)] + \
               [self.translate(reverse, start) for start in range(3)]

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
