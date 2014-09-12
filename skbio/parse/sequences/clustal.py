# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function

from skbio.io import RecordError
from skbio.parse.record import DelimitedSplitter


def _label_line_parser(record, splitter, strict=True):
    """Returns dict mapping list of data to labels, plus list with field order.

    Field order contains labels in order encountered in file.

    NOTE: doesn't care if lines are out of order in different blocks. This
    should never happen anyway, but it's possible that this behavior should
    be changed to tighten up validation.
    """
    labels = []
    result = {}
    for line in record:
        try:
            key, val = splitter(line.rstrip())
        except:
            if strict:
                raise RecordError(
                    "Failed to extract key and value from line %s" %
                    line)
            else:
                continue  # just skip the line if not strict

        if key in result:
            result[key].append(val)
        else:
            result[key] = [val]
            labels.append(key)
    return result, labels


def _is_clustal_seq_line(line):
    """Returns True if line starts with a non-blank character but not 'CLUSTAL'

    Useful for filtering other lines out of the file.
    """
    return line and (not line[0].isspace()) and\
        (not line.startswith('CLUSTAL')) and (not line.startswith('MUSCLE'))

last_space = DelimitedSplitter(None, -1)


def _delete_trailing_number(line):
    """Deletes trailing number from a line.

    WARNING: does not preserve internal whitespace when a number is removed!
    (converts each whitespace run to a single space). Returns the original
    line if it didn't end in a number.
    """
    pieces = line.split()
    try:
        int(pieces[-1])
        return ' '.join(pieces[:-1])
    except ValueError:  # no trailing numbers
        return line


def write_clustal(records, fh):
    clen = 60
    records = list(records)
    names, seqs = zip(*records)
    nameLen = max(map(len, names))
    seqLen = max(map(len, seqs))
    fh.write('CLUSTAL\n\n')
    for i in range(0, seqLen, clen):
        for label, seq in records:
            name = ('{:<%d}' % (nameLen)).format(label)
            fh.write("%s\t%s\t\n" % (name, seq[i:i+clen]))
        fh.write("\n")


def parse_clustal(record, strict=True):
    records = map(_delete_trailing_number,
                  filter(_is_clustal_seq_line, record))
    data, labels = _label_line_parser(records, last_space, strict)

    for key in labels:
        yield key, ''.join(data[key])
