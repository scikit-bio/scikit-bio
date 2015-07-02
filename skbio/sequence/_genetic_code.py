# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np
from future.builtins import range

from skbio.util._decorator import classproperty
from skbio._base import SkbioObject
from skbio.sequence import Protein, RNA
from skbio.sequence._base import ElasticLines


class GeneticCode(SkbioObject):
    _radix_multiplier = np.asarray([16, 4, 1], dtype=np.uint8)
    _byte_to_offset_map = {
        ord(b'U'): 0,
        ord(b'C'): 1,
        ord(b'A'): 2,
        ord(b'G'): 3
    }
    _start_stop_options = ['ignore', 'optional', 'require']

    @classmethod
    def from_ncbi(cls, table_id):
        if table_id not in _ncbi_genetic_codes:
            raise ValueError(
                "`table_id` must be one of %r, not %r"
                % (sorted(_ncbi_genetic_codes), table_id))
        return _ncbi_genetic_codes[table_id]

    @classproperty
    def reading_frames(cls):
        return [1, 2, 3, -1, -2, -3]

    @property
    def name(self):
        return self._name

    def __init__(self, amino_acids, starts, name=''):
        self._set_amino_acids(amino_acids)
        self._set_starts(starts)
        self._name = name

    def _set_amino_acids(self, amino_acids):
        amino_acids = Protein(amino_acids)

        if len(amino_acids) != 64:
            raise ValueError("`amino_acids` must be length 64, not %d"
                             % len(amino_acids))
        indices = (amino_acids.values == b'M').nonzero()[0]
        if indices.size < 1:
            raise ValueError("`amino_acids` must contain at least one M "
                             "(methionine) character")
        self._amino_acids = amino_acids
        self._m_character_codon = self._index_to_codon(indices[0])

    def _set_starts(self, starts):
        starts = Protein(starts)

        if len(starts) != 64:
            raise ValueError("`starts` must be length 64, not %d"
                             % len(starts))
        if ((starts.values == b'M').sum() + (starts.values == b'-').sum() !=
            len(starts)):
            # to prevent the user from accidentally swapping `starts` and
            # `amino_acids` and getting a translation back
            raise ValueError("`starts` may only contain M and - characters")

        self._starts = starts

        indices = (self._starts.values == b'M').nonzero()[0]
        codons = np.empty((indices.size, 3), dtype=np.uint8)
        for i, index in enumerate(indices):
            codons[i] = self._index_to_codon(index)
        self._start_codons = codons

    def _index_to_codon(self, index):
        codon = np.empty(3, dtype=np.uint8)
        for i, multiplier in enumerate(self._radix_multiplier):
            offset, index = divmod(index, multiplier)
            codon[i] = offset
        return codon

    def __str__(self):
        return self._build_repr(include_name=False)

    def __repr__(self):
        return self._build_repr(include_name=True)

    def _build_repr(self, include_name):
        lines = ElasticLines()

        if include_name:
            name_line = self.__class__.__name__
            if len(self.name) > 0:
                name_line += ' (%s)' % self.name
            lines.add_line(name_line)
            lines.add_separator()

        lines.add_line('  AAs  = %s' % str(self._amino_acids))
        lines.add_line('Starts = %s' % str(self._starts))
        base1 = 'U' * 16 + 'C' * 16 + 'A' * 16 + 'G' * 16
        lines.add_line('Base1  = %s' % base1)
        base2 = ('U' * 4 + 'C' * 4 + 'A' * 4 + 'G' * 4) * 4
        lines.add_line('Base2  = %s' % base2)
        base3 = 'UCAG' * 16
        lines.add_line('Base3  = %s' % base3)

        return lines.to_str()

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        if str(self._amino_acids) != str(other._amino_acids):
            return False
        if str(self._starts) != str(other._starts):
            return False
        return True

    def __ne__(self, other):
        return not (self == other)

    def translate(self, sequence, reading_frame=1, start='ignore',
                  stop='ignore'):
        """

        Parameters
        ----------
        sequence : RNA
            RNA sequence to translate.
        reading_frame : {1, 2, 3, -1, -2, -3}
            Reading frame to use. The number indicates the base to start
            translation on. If negative, will perform a reverse complement
            first.
        start : {'ignore', 'optional', 'require'}
            If ``True``, translation will begin at the first start codon in the
            reading frame, ignoring all bases prior. 'M' will always be the
            first amino acid in the translated sequence, even if the original
            start codon coded for a different amino acid (this behavior most
            closely matches what happens biologically). If ``False``,
            translation will start from the position indicated by the reading
            frame, regardless of the presence of a start codon.
        stop : {'ignore', 'optional', 'require'}
            If ``True``, translation will terminate at the first stop codon.
            The stop codon will not be included in the translated sequence. If
            ``False``, translation will terminate at the last codon in the
            sequence, even if it is not a stop codon.

        Returns
        -------
        Protein
            Translated protein sequence.

        """
        if not isinstance(sequence, RNA):
            raise TypeError("Sequence to translate must be RNA, not %s" %
                            type(sequence).__name__)

        if reading_frame not in self.reading_frames:
            raise ValueError("`reading_frame` must be one of %r, not %r" %
                             (self.reading_frames, reading_frame))

        for name, value in ('start', start), ('stop', stop):
            if value not in self._start_stop_options:
                raise ValueError("`%s` must be one of %r, not %r" %
                                 (name, self._start_stop_options, value))

        if sequence.has_gaps():
            raise ValueError("scikit-bio does not support translation of "
                             "gapped sequences.")

        if sequence.has_degenerates():
            raise NotImplementedError("scikit-bio does not currently support "
                                      "translation of degenerate sequences.")

        offset = abs(reading_frame) - 1
        if reading_frame < 0:
            sequence = sequence.reverse_complement()

        data = sequence.values[offset:].view(np.uint8).copy()
        data = data[:data.size // 3 * 3].reshape((-1, 3))

        for byte in self._byte_to_offset_map:
            data[data == byte] = self._byte_to_offset_map[byte]

        if start in {'require', 'optional'}:
            start_codon_index = data.shape[0]
            for start_codon in self._start_codons:
                indices = np.all(data == start_codon, axis=1).nonzero()[0]

                if indices.size > 0:
                    first_index = indices[0]
                    if first_index < start_codon_index:
                        start_codon_index = first_index

            if start_codon_index != data.shape[0]:
                data = data[start_codon_index:]
                # replace start codon with codon coding for M
                data[0] = self._m_character_codon
            elif start == 'require':
                self._raise_require_error('start', reading_frame)

        index = (data * self._radix_multiplier).sum(axis=1)
        translated = self._amino_acids.values[index]

        if stop in {'require', 'optional'}:
            stop_codon_indices = (translated == b'*').nonzero()[0]
            if stop_codon_indices.size > 0:
                translated = translated[:stop_codon_indices[0]]
            elif stop == 'require':
                self._raise_require_error('stop', reading_frame)

        metadata = None
        if sequence.has_metadata():
            metadata = sequence.metadata

        return Protein(translated, metadata=metadata)

    def _raise_require_error(self, name, reading_frame):
        raise ValueError(
            "Sequence does not contain a %s codon in the "
            "current reading frame (`reading_frame=%d`). Presence "
            "of a %s codon is required with `%s='require'`"
            % (name, reading_frame, name, name))

    def translate_six_frames(self, sequence, start='ignore', stop='ignore'):
        rc = sequence.reverse_complement()

        for reading_frame in range(1, 4):
            yield self.translate(sequence, reading_frame=reading_frame,
                                 start=start, stop=stop)
        for reading_frame in range(1, 4):
            yield self.translate(rc, reading_frame=reading_frame,
                                 start=start, stop=stop)


_ncbi_genetic_codes = {
    1: GeneticCode(
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '---M---------------M---------------M----------------------------',
        'Standard'),
    2: GeneticCode(
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG',
        '--------------------------------MMMM---------------M------------',
        'Vertebrate Mitochondrial'),
    3: GeneticCode(
        'FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '----------------------------------MM----------------------------',
        'Yeast Mitochondrial'),
    4: GeneticCode(
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '--MM---------------M------------MMMM---------------M------------',
        'Mold, Protozoan, and Coelenterate Mitochondrial, and '
        'Mycoplasma/Spiroplasma'),
    5: GeneticCode(
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG',
        '---M----------------------------MMMM---------------M------------',
        'Invertebrate Mitochondrial'),
    6: GeneticCode(
        'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '-----------------------------------M----------------------------',
        'Ciliate, Dasycladacean and Hexamita Nuclear'),
    9: GeneticCode(
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        '-----------------------------------M---------------M------------',
        'Echinoderm and Flatworm Mitochondrial'),
    10: GeneticCode(
        'FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '-----------------------------------M----------------------------',
        'Euplotid Nuclear'),
    11: GeneticCode(
        'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '---M---------------M------------MMMM---------------M------------',
        'Bacterial, Archaeal and Plant Plastid'),
    12: GeneticCode(
        'FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '-------------------M---------------M----------------------------',
        'Alternative Yeast Nuclear'),
    13: GeneticCode(
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG',
        '---M------------------------------MM---------------M------------',
        'Ascidian Mitochondrial'),
    14: GeneticCode(
        'FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        '-----------------------------------M----------------------------',
        'Alternative Flatworm Mitochondrial'),
    16: GeneticCode(
        'FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '-----------------------------------M----------------------------',
        'Chlorophycean Mitochondrial'),
    21: GeneticCode(
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG',
        '-----------------------------------M---------------M------------',
        'Trematode Mitochondrial'),
    22: GeneticCode(
        'FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '-----------------------------------M----------------------------',
        'Scenedesmus obliquus Mitochondrial'),
    23: GeneticCode(
        'FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '--------------------------------M--M---------------M------------',
        'Thraustochytrium Mitochondrial'),
    24: GeneticCode(
        'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG',
        '---M---------------M---------------M---------------M------------',
        'Pterobranchia Mitochondrial'),
    25: GeneticCode(
        'FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        '---M-------------------------------M---------------M------------',
        'Candidate Division SR1 and Gracilibacteria')
}
