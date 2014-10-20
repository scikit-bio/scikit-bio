# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.standard_library import hooks
with hooks():
    from itertools import zip_longest

from skbio.alignment import SequenceCollection, Alignment
from skbio.io import (register_reader, register_writer, register_sniffer,
                      FASTAQUALFormatError)
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)

from ._base import (
    _fasta_or_qual_sniffer, _fasta_or_qual_to_generator,
    _fasta_or_fasta_qual_to_sequence, _generator_to_fasta_or_fasta_qual)


@register_sniffer('qual')
def _qual_sniffer(fh):
    return _fasta_or_qual_sniffer(fh, 'qual')


# TODO is this is correct way to handle the constructor kwarg in a compound
# format? should there be a check that the two values are exactly the same?
@register_reader(['fasta', 'qual'])
def _fasta_qual_to_generator(fasta_fh, qual_fh,
                             constructor=(BiologicalSequence,
                                          BiologicalSequence)):
    # TODO is this try/finally necessary anymore?
    try:
        fasta_gen = _fasta_or_qual_to_generator(fasta_fh, 'fasta')
        qual_gen = _fasta_or_qual_to_generator(qual_fh, 'qual')

        for fasta_rec, qual_rec in zip_longest(fasta_gen, qual_gen,
                                               fillvalue=None):
            if fasta_rec is None:
                raise FASTAQUALFormatError(
                    "QUAL file has more records than FASTA file.")
            if qual_rec is None:
                raise FASTAQUALFormatError(
                    "FASTA file has more records than QUAL file.")

            fasta_seq, fasta_id, fasta_desc = fasta_rec
            qual_scores, qual_id, qual_desc = qual_rec

            if fasta_id != qual_id:
                raise FASTAQUALFormatError(
                    "IDs do not match between FASTA and QUAL records: %r != %r"
                    % (fasta_id, qual_id))
            if fasta_desc != qual_desc:
                raise FASTAQUALFormatError(
                    "Descriptions do not match between FASTA and QUAL "
                    "records: %r != %r" % (fasta_desc, qual_desc))

            yield constructor[0](fasta_seq, id=fasta_id,
                                 description=fasta_desc, quality=qual_scores)
    finally:
        fasta_gen.close()
        qual_gen.close()


@register_reader(['fasta', 'qual'], BiologicalSequence)
def _fasta_qual_to_biological_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0],
        (BiologicalSequence, BiologicalSequence), _fasta_qual_to_generator,
        FASTAQUALFormatError, 'FASTA/QUAL')


@register_reader(['fasta', 'qual'], NucleotideSequence)
def _fasta_qual_to_nucleotide_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0],
        (NucleotideSequence, NucleotideSequence), _fasta_qual_to_generator,
        FASTAQUALFormatError, 'FASTA/QUAL')


@register_reader(['fasta', 'qual'], DNASequence)
def _fasta_qual_to_dna_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0], (DNASequence, DNASequence),
        _fasta_qual_to_generator, FASTAQUALFormatError, 'FASTA/QUAL')


@register_reader(['fasta', 'qual'], RNASequence)
def _fasta_qual_to_rna_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0], (RNASequence, RNASequence),
        _fasta_qual_to_generator, FASTAQUALFormatError, 'FASTA/QUAL')


@register_reader(['fasta', 'qual'], ProteinSequence)
def _fasta_qual_to_protein_sequence(fasta_fh, qual_fh, seq_num=(1, 1)):
    return _fasta_or_fasta_qual_to_sequence(
        (fasta_fh, qual_fh), seq_num[0], (ProteinSequence, ProteinSequence),
        _fasta_qual_to_generator, FASTAQUALFormatError, 'FASTA/QUAL')


@register_reader(['fasta', 'qual'], SequenceCollection)
def _fasta_qual_to_sequence_collection(fasta_fh, qual_fh,
                                       constructor=(BiologicalSequence,
                                                    BiologicalSequence)):
    return SequenceCollection(
        list(_fasta_qual_to_generator(fasta_fh, qual_fh,
                                      constructor=constructor)))


@register_reader(['fasta', 'qual'], Alignment)
def _fasta_qual_to_alignment(fasta_fh, qual_fh,
                             constructor=(BiologicalSequence,
                                          BiologicalSequence)):
    return Alignment(
        list(_fasta_qual_to_generator(fasta_fh, qual_fh,
                                      constructor=constructor)))


@register_writer(['fasta', 'qual'])
def _generator_to_fasta_qual(obj, fasta_fh, qual_fh,
                             id_whitespace_replacement=('_', '_'),
                             description_newline_replacement=(' ', ' '),
                             max_width=(None, None)):
    _generator_to_fasta_or_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement[0],
        description_newline_replacement[0], max_width[0])


@register_writer(['fasta', 'qual'], BiologicalSequence)
def _biological_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                       id_whitespace_replacement=('_', '_'),
                                       description_newline_replacement=(' ',
                                                                        ' '),
                                       max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], NucleotideSequence)
def _nucleotide_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                       id_whitespace_replacement=('_', '_'),
                                       description_newline_replacement=(' ',
                                                                        ' '),
                                       max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], DNASequence)
def _dna_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                id_whitespace_replacement=('_', '_'),
                                description_newline_replacement=(' ', ' '),
                                max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], RNASequence)
def _rna_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                id_whitespace_replacement=('_', '_'),
                                description_newline_replacement=(' ', ' '),
                                max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], ProteinSequence)
def _protein_sequence_to_fasta_qual(obj, fasta_fh, qual_fh,
                                    id_whitespace_replacement=('_', '_'),
                                    description_newline_replacement=(' ', ' '),
                                    max_width=(None, None)):
    _sequence_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], SequenceCollection)
def _sequence_collection_to_fasta_qual(obj, fasta_fh, qual_fh,
                                       id_whitespace_replacement=('_', '_'),
                                       description_newline_replacement=(' ',
                                                                        ' '),
                                       max_width=(None, None)):
    _sequences_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


@register_writer(['fasta', 'qual'], Alignment)
def _alignment_to_fasta_qual(obj, fasta_fh, qual_fh,
                             id_whitespace_replacement=('_', '_'),
                             description_newline_replacement=(' ', ' '),
                             max_width=(None, None)):
    _sequences_to_fasta_qual(
        obj, fasta_fh, qual_fh, id_whitespace_replacement,
        description_newline_replacement, max_width)


def _sequence_to_fasta_qual(obj, fasta_fh, qual_fh, id_whitespace_replacement,
                            description_newline_replacement, max_width):
    def seq_gen():
        yield obj

    _generator_to_fasta_qual(
        seq_gen(), fasta_fh, qual_fh,
        id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width)


def _sequences_to_fasta_qual(obj, fasta_fh, qual_fh, id_whitespace_replacement,
                             description_newline_replacement, max_width):
    def seq_gen():
        for seq in obj:
            yield seq

    _generator_to_fasta_qual(
        seq_gen(), fasta_fh, qual_fh,
        id_whitespace_replacement=id_whitespace_replacement,
        description_newline_replacement=description_newline_replacement,
        max_width=max_width)
