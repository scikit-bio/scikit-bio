# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.alignment import SequenceCollection, Alignment
from skbio.io import (register_reader, register_writer, register_sniffer,
                      FASTAFormatError)
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)

from ._base import (
    _fasta_or_qual_sniffer, _fasta_or_qual_to_generator,
    _fasta_or_fasta_qual_to_sequence, _generator_to_fasta_or_fasta_qual)


@register_sniffer('fasta')
def _fasta_sniffer(fh):
    # Strategy:
    #   If the file appears to be valid FASTA, check if it's also valid QUAL.
    #   If it is, do not identify as FASTA since the QUAL sniffer is more
    #   specific (i.e., it attempts to parse the quality scores as integers).
    valid_fasta = _fasta_or_qual_sniffer(fh, 'fasta')[0]
    if valid_fasta:
        fh.seek(0)
        valid_fasta = not _fasta_or_qual_sniffer(fh, 'qual')[0]
    return valid_fasta, {}


@register_reader('fasta')
def _fasta_to_generator(fh, constructor=BiologicalSequence):
    for seq, id_, desc in _fasta_or_qual_to_generator(fh, format='fasta'):
        yield constructor(seq, id=id_, description=desc)


@register_reader('fasta', BiologicalSequence)
def _fasta_to_biological_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence(
        (fh,), seq_num, BiologicalSequence, _fasta_to_generator,
        FASTAFormatError, 'FASTA')


@register_reader('fasta', NucleotideSequence)
def _fasta_to_nucleotide_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence(
        (fh,), seq_num, NucleotideSequence, _fasta_to_generator,
        FASTAFormatError, 'FASTA')


@register_reader('fasta', DNASequence)
def _fasta_to_dna_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence(
        (fh,), seq_num, DNASequence, _fasta_to_generator, FASTAFormatError,
        'FASTA')


@register_reader('fasta', RNASequence)
def _fasta_to_rna_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence(
        (fh,), seq_num, RNASequence, _fasta_to_generator, FASTAFormatError,
        'FASTA')


@register_reader('fasta', ProteinSequence)
def _fasta_to_protein_sequence(fh, seq_num=1):
    return _fasta_or_fasta_qual_to_sequence(
        (fh,), seq_num, ProteinSequence, _fasta_to_generator, FASTAFormatError,
        'FASTA')


@register_reader('fasta', SequenceCollection)
def _fasta_to_sequence_collection(fh, constructor=BiologicalSequence):
    return SequenceCollection(
        list(_fasta_to_generator(fh, constructor=constructor)))


@register_reader('fasta', Alignment)
def _fasta_to_alignment(fh, constructor=BiologicalSequence):
    return Alignment(
        list(_fasta_to_generator(fh, constructor=constructor)))


@register_writer('fasta')
def _generator_to_fasta(obj, fh, id_whitespace_replacement='_',
                        description_newline_replacement=' ', max_width=None):
    _generator_to_fasta_or_fasta_qual(
        obj, fh, None, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer('fasta', BiologicalSequence)
def _biological_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', NucleotideSequence)
def _nucleotide_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', DNASequence)
def _dna_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                           description_newline_replacement=' ',
                           max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', RNASequence)
def _rna_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                           description_newline_replacement=' ',
                           max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', ProteinSequence)
def _protein_sequence_to_fasta(obj, fh, id_whitespace_replacement='_',
                               description_newline_replacement=' ',
                               max_width=None):
    _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width)


@register_writer('fasta', SequenceCollection)
def _sequence_collection_to_fasta(obj, fh, id_whitespace_replacement='_',
                                  description_newline_replacement=' ',
                                  max_width=None):
    _sequences_to_fasta(obj, fh, id_whitespace_replacement,
                        description_newline_replacement, max_width)


@register_writer('fasta', Alignment)
def _alignment_to_fasta(obj, fh, id_whitespace_replacement='_',
                        description_newline_replacement=' ', max_width=None):
    _sequences_to_fasta(obj, fh, id_whitespace_replacement,
                        description_newline_replacement, max_width)


def _sequence_to_fasta(obj, fh, id_whitespace_replacement,
                       description_newline_replacement, max_width):
    def seq_gen():
        yield obj

    _generator_to_fasta(
        seq_gen(), fh, id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width)


def _sequences_to_fasta(obj, fh, id_whitespace_replacement,
                        description_newline_replacement, max_width):
    def seq_gen():
        for seq in obj:
            yield seq

    _generator_to_fasta(
        seq_gen(), fh, id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width)
