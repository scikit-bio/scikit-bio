# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.builtins import range, zip

import re
import textwrap

import numpy as np

from skbio.io import FASTAFormatError, QUALFormatError, FASTAQUALFormatError
from skbio.io._base import _chunk_str
from skbio.util import cardinal_to_ordinal


def _fasta_or_qual_sniffer(fh, format):
    # Strategy:
    #   Read up to 10 records. If at least one record is read (i.e. the file
    #   isn't empty) and no errors are thrown during reading, assume the file
    #   is in FASTA or QUAL format.

    # TODO is the finally block still necessary?
    try:
        not_empty = False
        gen = _fasta_or_qual_to_generator(fh, format=format)
        for _ in zip(range(10), gen):
            not_empty = True
        return not_empty, {}
    except (FASTAFormatError, QUALFormatError):
        return False, {}
    finally:
        gen.close()


def _fasta_or_qual_to_generator(fh, format):
    if format == 'fasta':
        data_parser = _parse_sequence_data
        error_type = FASTAFormatError
        format_label = 'FASTA'
    else:
        data_parser = _parse_quality_scores
        error_type = QUALFormatError
        format_label = 'QUAL'

    line = next(fh)
    # header check inlined here and below for performance
    if line.startswith('>'):
        id_, desc = _parse_header(line)
    else:
        raise error_type(
            "Found line without a %s header:\n%s" % (format_label, line))

    data_chunks = []
    for line in fh:
        if line.startswith('>'):
            # new header, so yield current record and reset state
            yield data_parser(data_chunks), id_, desc
            data_chunks = []
            id_, desc = _parse_header(line)
        else:
            line = line.strip()
            if line:
                data_chunks.append(line)
            else:
                raise error_type("Found blank or whitespace-only line in "
                                 "%s-formatted file." % format_label)
    # yield last record in file
    yield data_parser(data_chunks), id_, desc


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


def _parse_sequence_data(chunks):
    if not chunks:
        raise FASTAFormatError("Found FASTA header without sequence data.")
    return ''.join(chunks)


def _parse_quality_scores(chunks):
    if not chunks:
        raise QUALFormatError("Found QUAL header without quality scores.")

    qual_str = ' '.join(chunks)
    try:
        return np.asarray(qual_str.split(), dtype=int)
    except ValueError:
        raise QUALFormatError(
            "Could not convert quality scores to integers:\n%s" % qual_str)


def _fasta_or_fasta_qual_to_sequence(fh, seq_num, constructor, seq_generator,
                                     error_type, format_label):
    if seq_num < 1:
        raise error_type(
            "Invalid sequence number (seq_num=%d). seq_num must be between 1 "
            "and the number of sequences in the %s-formatted file (inclusive)."
            % (seq_num, format_label))

    seq_idx = seq_num - 1
    seq = None
    try:
        gen = seq_generator(*fh, constructor=constructor)
        for idx, curr_seq in enumerate(gen):
            if idx == seq_idx:
                seq = curr_seq
                break
    finally:
        gen.close()

    if seq is None:
        raise error_type(
            "Reached end of %s-formatted file before finding %s biological "
            "sequence." % (format_label, cardinal_to_ordinal(seq_num)))
    return seq


def _generator_to_fasta_or_fasta_qual(obj, fasta_fh, qual_fh,
                                      id_whitespace_replacement,
                                      description_newline_replacement,
                                      max_width):
    if qual_fh is None:
        error_type = FASTAFormatError
        format_label = 'FASTA'
    else:
        error_type = FASTAQUALFormatError
        format_label = 'FASTA/QUAL'
        if max_width is not None:
            # define text wrapper for quality scores here for efficiency.
            # textwrap docs recommend reusing a TextWrapper instance when it is
            # used many times. configure text wrapper to never break words
            # (i.e., integer quality scores)
            qual_wrapper = textwrap.TextWrapper(width=max_width,
                                                break_long_words=False,
                                                break_on_hyphens=False)

    if ((id_whitespace_replacement is not None and
         '\n' in id_whitespace_replacement) or
        (description_newline_replacement is not None and
         '\n' in description_newline_replacement)):
        raise error_type(
            "Newline character (\\n) cannot be used to replace whitespace in "
            "biological sequence IDs, nor to replace newlines in biological "
            "sequence descriptions. Otherwise, the %s-formatted file will be "
            "invalid." % format_label)
    ws_pattern = re.compile(r'\s')
    nl_pattern = re.compile(r'\n')

    for idx, seq in enumerate(obj):
        if len(seq) < 1:
            raise error_type(
                "Cannot write %s biological sequence in %s format because it "
                "does not contain any characters (i.e., it is an empty/blank "
                "sequence). Empty sequences are not supported in the %s file "
                "format." % (cardinal_to_ordinal(idx + 1), format_label,
                             format_label))

        id_ = seq.id
        if id_whitespace_replacement is not None:
            id_ = re.sub(ws_pattern, id_whitespace_replacement, id_)

        desc = seq.description
        if description_newline_replacement is not None:
            desc = re.sub(nl_pattern, description_newline_replacement, desc)

        if desc:
            header = '%s %s' % (id_, desc)
        else:
            header = id_

        seq_str = str(seq)
        if max_width is not None:
            seq_str = _chunk_str(seq_str, max_width, '\n')

        fasta_fh.write('>%s\n%s\n' % (header, seq_str))

        if qual_fh is not None:
            if not seq.has_quality():
                raise error_type(
                    "Cannot write %s biological sequence in %s format because "
                    "it does not have quality scores associated with it." %
                    (cardinal_to_ordinal(idx + 1), format_label))

            qual_str = ' '.join(np.asarray(seq.quality, dtype=np.str))
            if max_width is not None:
                qual_str = qual_wrapper.fill(qual_str)

            qual_fh.write('>%s\n%s\n' % (header, qual_str))
