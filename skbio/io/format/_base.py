# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
import warnings

import numpy as np

from skbio.util import cardinal_to_ordinal

_whitespace_regex = re.compile(r"\s")
_newline_regex = re.compile(r"\n")


def _decode_qual_to_phred(qual_str, variant=None, phred_offset=None):
    phred_offset, phred_range = _get_phred_offset_and_range(
        variant,
        phred_offset,
        [
            "Must provide either `variant` or `phred_offset` in order to decode "
            "quality scores.",
            "Decoding Solexa quality scores is not currently supported, "
            "as quality scores are always stored as Phred scores in "
            "scikit-bio. Please see the following scikit-bio issue to "
            "track progress on this:\n\t"
            "https://github.com/scikit-bio/scikit-bio/issues/719",
        ],
    )
    qual = np.frombuffer(qual_str.encode("ascii"), dtype=np.uint8) - phred_offset

    if np.any((qual > phred_range[1]) | (qual < phred_range[0])):
        raise ValueError(
            "Decoded Phred score is out of range [%d, %d]."
            % (phred_range[0], phred_range[1])
        )

    return qual


def _encode_phred_to_qual(phred, variant=None, phred_offset=None):
    phred_offset, phred_range = _get_phred_offset_and_range(
        variant,
        phred_offset,
        [
            "Must provide either `variant` or `phred_offset` in order to encode "
            "Phred scores.",
            "Encoding Solexa quality scores is not currently supported. "
            "Please see the following scikit-bio issue to track progress "
            "on this:\n\t"
            "https://github.com/scikit-bio/scikit-bio/issues/719",
        ],
    )

    qual_chars = []
    for score in phred:
        if score < phred_range[0]:
            raise ValueError(
                "Phred score %d is out of range [%d, %d]."
                % (score, phred_range[0], phred_range[1])
            )
        if score > phred_range[1]:
            warnings.warn(
                "Phred score %d is out of targeted range [%d, %d]. Converting "
                "to %d." % (score, phred_range[0], phred_range[1], phred_range[1]),
                UserWarning,
            )
            score = phred_range[1]
        qual_chars.append(chr(score + phred_offset))
    return "".join(qual_chars)


def _get_phred_offset_and_range(variant, phred_offset, errors):
    if variant is None and phred_offset is None:
        raise ValueError(errors[0])
    if variant is not None and phred_offset is not None:
        raise ValueError("Cannot provide both `variant` and `phred_offset`.")

    if variant is not None:
        if variant == "sanger":
            phred_offset = 33
            phred_range = (0, 93)
        elif variant == "illumina1.3":
            phred_offset = 64
            phred_range = (0, 62)
        elif variant == "illumina1.8":
            phred_offset = 33
            phred_range = (0, 62)
        elif variant == "solexa":
            phred_offset = 64
            phred_range = (-5, 62)
            raise ValueError(errors[1])
        else:
            raise ValueError("Unrecognized variant %r." % variant)
    else:
        if not (33 <= phred_offset <= 126):
            raise ValueError(
                "`phred_offset` %d is out of printable ASCII character range."
                % phred_offset
            )
        phred_range = (0, 126 - phred_offset)

    return phred_offset, phred_range


def _get_nth_sequence(generator, seq_num):
    # i is set to None so that an empty generator will not result in an
    # undefined variable when compared to seq_num.
    i = None
    if seq_num is None or seq_num < 1:
        raise ValueError(
            "Invalid sequence number (`seq_num`=%s). `seq_num`"
            " must be between 1 and the number of sequences in"
            " the file." % str(seq_num)
        )
    try:
        for i, seq in zip(range(1, seq_num + 1), generator):
            pass
    finally:
        generator.close()

    if i == seq_num:
        return seq
    raise ValueError(
        "Reached end of file before finding the %s sequence."
        % cardinal_to_ordinal(seq_num)
    )


def _parse_fasta_like_header(line):
    id_ = ""
    desc = ""
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


def _format_fasta_like_records(
    generator,
    id_whitespace_replacement,
    description_newline_replacement,
    require_qual,
    lowercase=None,
):
    if (
        id_whitespace_replacement is not None and "\n" in id_whitespace_replacement
    ) or (
        description_newline_replacement is not None
        and "\n" in description_newline_replacement
    ):
        raise ValueError(
            "Newline character (\\n) cannot be used to replace whitespace in "
            "sequence IDs, nor to replace newlines in sequence descriptions."
        )

    for idx, seq in enumerate(generator):
        if len(seq) < 1:
            raise ValueError(
                "%s sequence does not contain any characters (i.e., it is an "
                "empty/blank sequence). Writing empty sequences is not "
                "supported." % cardinal_to_ordinal(idx + 1)
            )

        if "id" in seq.metadata:
            id_ = "%s" % seq.metadata["id"]
        else:
            id_ = ""

        if id_whitespace_replacement is not None:
            id_ = _whitespace_regex.sub(id_whitespace_replacement, id_)

        if "description" in seq.metadata:
            desc = "%s" % seq.metadata["description"]
        else:
            desc = ""

        if description_newline_replacement is not None:
            desc = _newline_regex.sub(description_newline_replacement, desc)

        if desc:
            header = "%s %s" % (id_, desc)
        else:
            header = id_

        if require_qual and "quality" not in seq.positional_metadata:
            raise ValueError(
                "Cannot write %s sequence because it does not have quality "
                "scores associated with it." % cardinal_to_ordinal(idx + 1)
            )

        qual = None
        if "quality" in seq.positional_metadata:
            qual = seq.positional_metadata["quality"].values

        if lowercase is not None:
            seq_str = seq.lowercase(lowercase)
        else:
            seq_str = str(seq)
        yield header, "%s" % seq_str, qual


def _line_generator(fh, skip_blanks=False, strip=True):
    for line in fh:
        if strip:
            line = line.strip()
        skip = False
        if skip_blanks:
            skip = line.isspace() or not line
        if not skip:
            yield line


def _too_many_blanks(fh, max_blanks):
    count = 0
    too_many = False
    for line in _line_generator(fh, skip_blanks=False):
        if line:
            break
        else:
            count += 1
            if count > max_blanks:
                too_many = True
                break
    fh.seek(0)
    return too_many
