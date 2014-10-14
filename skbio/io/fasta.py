"""
FASTA format (:mod:`skbio.io.fasta`)
====================================

.. currentmodule:: skbio.io.fasta

TODO add description

Format Support
--------------
**Has Sniffer: Yes**

+----------+----------+------------------------------------------------------+
|**Reader**|**Writer**|                   **Object Class**                   |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |generator of ``BiologicalSequence`` objects           |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.alignment.SequenceCollection`             |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.alignment.Alignment`                      |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.BiologicalSequence`              |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.NucleotideSequence`              |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.DNASequence`                     |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.RNASequence`                     |
+----------+----------+------------------------------------------------------+
|Yes       |Yes       |:mod:`skbio.sequence.ProteinSequence`                 |
+----------+----------+------------------------------------------------------+

Format Specification
--------------------
TODO add format specification

Caveats
^^^^^^^

Format Parameters
-----------------
TODO add format parameters

Examples
--------
TODO add examples

References
----------
.. [1] Lipman, DJ; Pearson, WR (1985). "Rapid and sensitive protein similarity
   searches". Science 227 (4693): 1435–41.
.. [2] http://en.wikipedia.org/wiki/FASTA_format
.. [3] http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import re

from skbio.io import (register_reader, register_writer, register_sniffer,
                      FASTAFormatError)
from skbio.io._base import _chunk_str
from skbio.alignment import SequenceCollection, Alignment
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)
from skbio.util import cardinal_to_ordinal


@register_sniffer('fasta')
def _fasta_sniffer(fh):
    # Strategy:
    #   Read up to 10 sequences (records). If at least one sequence is read
    #   (i.e. the file isn't empty) and no errors are thrown during reading,
    #   assume the file is in FASTA format.
    try:
        gen = _fasta_to_generator(fh)
        num_seqs = 0
        for _ in range(10):
            try:
                next(gen)
            except FASTAFormatError:
                return False, {}
            except StopIteration:
                break
            num_seqs += 1

        if num_seqs < 1:
            return False, {}
        else:
            return True, {}
    finally:
        gen.close()


@register_reader('fasta')
def _fasta_to_generator(fh, constructor=BiologicalSequence):
    line = next(fh)
    if _is_header(line):
        id_, desc = _parse_header(line)
    else:
        raise FASTAFormatError("Found line without a FASTA header:\n%s" % line)

    seq_chunks = []
    for line in fh:
        if _is_header(line):
            # new header, so yield current sequence and reset state
            yield _construct_sequence(constructor, seq_chunks, id_, desc)
            seq_chunks = []
            id_, desc = _parse_header(line)
        else:
            line = line.strip()
            if line:
                seq_chunks.append(line)
            else:
                raise FASTAFormatError("Found blank or whitespace-only line "
                                       "in FASTA-formatted file.")

    # yield last sequence in file
    yield _construct_sequence(constructor, seq_chunks, id_, desc)


@register_reader('fasta', BiologicalSequence)
def _fasta_to_biological_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, BiologicalSequence)


@register_reader('fasta', NucleotideSequence)
def _fasta_to_nucleotide_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, NucleotideSequence)


@register_reader('fasta', DNASequence)
def _fasta_to_dna_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, DNASequence)


@register_reader('fasta', RNASequence)
def _fasta_to_rna_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, RNASequence)


@register_reader('fasta', ProteinSequence)
def _fasta_to_protein_sequence(fh, seq_num=1):
    return _fasta_to_sequence(fh, seq_num, ProteinSequence)


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
    if '\n' in id_whitespace_replacement or \
            '\n' in description_newline_replacement:
        raise FASTAFormatError(
            "Newline character (\\n) cannot be used to replace whitespace in "
            "biological sequence IDs, nor to replace newlines in biological "
            "sequence descriptions. Otherwise, the FASTA-formatted file will "
            "be invalid.")
    ws_pattern = re.compile(r'\s')
    nl_pattern = re.compile(r'\n')

    for idx, seq in enumerate(obj):
        if len(seq) < 1:
            raise FASTAFormatError(
                "Cannot write %s biological sequence in FASTA format because "
                "it does not contain any characters (i.e., it is an "
                "empty/blank sequence). Empty sequences are not supported in "
                "the FASTA file format." % cardinal_to_ordinal(idx + 1))

        id_ = re.sub(ws_pattern, id_whitespace_replacement, seq.id)
        desc = re.sub(nl_pattern, description_newline_replacement,
                      seq.description)

        if desc:
            header = '%s %s' % (id_, desc)
        else:
            header = id_

        seq_str = str(seq)

        if max_width is not None:
            seq_str = _chunk_str(seq_str, max_width, '\n')

        fh.write('>%s\n%s\n' % (header, seq_str))


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


def _is_header(line):
    return line.startswith('>')


def _parse_header(line):
    id_ = ''
    desc = ''
    header = line[1:].rstrip()
    if header:
        if header[0].isspace():
            # no id
            desc = header.lstrip()
        else:
            header_tokens = header.split(None, 1)
            if len(header_tokens) == 1:
                # no description
                id_ = header_tokens[0]
            else:
                id_, desc = header_tokens
    return id_, desc


def _construct_sequence(constructor, seq_chunks, id_, description):
    if not seq_chunks:
        raise FASTAFormatError("Found FASTA header without sequence data.")
    return constructor(''.join(seq_chunks), id=id_, description=description)


def _fasta_to_sequence(fh, seq_num, constructor):
    if seq_num < 1:
        raise FASTAFormatError(
            "Invalid sequence number (seq_num=%d). seq_num must be between 1 "
            "and the number of sequences in the FASTA-formatted file "
            "(inclusive)." % seq_num)

    seq_idx = seq_num - 1
    seq = None
    gen = _fasta_to_generator(fh, constructor=constructor)
    for idx, curr_seq in enumerate(gen):
        if idx == seq_idx:
            seq = curr_seq
            break
    # TODO is this necessary? what if an error is raised within the generator?
    gen.close()

    if seq is None:
        raise FASTAFormatError(
            "Reached end of FASTA-formatted file before finding %s biological "
            "sequence." % cardinal_to_ordinal(seq_num))
    return seq


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
