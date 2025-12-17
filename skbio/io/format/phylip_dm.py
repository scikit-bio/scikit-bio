"""PHYLIP distance matrix format (:mod:`skbio.io.format.phylip_dm`)
================================================================

.. currentmodule:: skbio.io.format.phylip_dm

The PHYLIP file format can store a distance matrix.

An example PHYLIP-formatted distance matrix::

    4
    Seq1      0.0000 1.6866 1.7198 1.6606
    Seq2      1.6866 0.0000 1.5232 1.4841
    Seq3      1.7198 1.5232 0.0000 0.7115
    Seq4      1.6606 1.4841 0.7115 0.0000

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.stats.distance.DistanceMatrix`                     |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
PHYLIP format is a plain text format containing exactly two sections: a header
describing the number of sequences in the matrix, followed by the matrix itself.

The format described here is "relaxed" PHYLIP. Relaxed PHYLIP does not require
that each sequence identifier is exactly 10 characters long like strict PHYLIP does.
Sequence identifiers and matrix data are separated by whitespace (spaces or tabs).

The format described here is "sequential" format. The original PHYLIP specification
describes both sequential and interleaved formats.

.. note:: scikit-bio currently supports reading and writing relaxed, sequential
   PHYLIP-formatted files. Strict and/or interleaved PHYLIP formats are not
   supported.

Header Section
^^^^^^^^^^^^^^
The header consists of a single line with a single positive integer (``n``) that
specifies the number of sequences in the matrix. This **must** be the first line
in the file. The integer may be preceded by optional whitespace.

.. note:: scikit-bio writes the PHYLIP format header without preceding spaces.

   PHYLIP format does not support blank line(s) between the header and the matrix.

Matrix Section
^^^^^^^^^^^^^^
The matrix section immediately follows the header. It consists of ``n`` lines (rows),
one for each sequence. Each row consists of a sequence identifier (ID) followed by
the distance values for that sequence, separated by whitespace.

**Square matrices**: Each row contains ``n`` distance values (the full row of the
distance matrix).

**Lower triangular matrices**: Row ``i`` contains ``i`` distance values (only the
values below the diagonal). The first row contains no distance values, just the ID.

Sequence IDs
^^^^^^^^^^^^
- IDs can have arbitrary length (relaxed format)
- IDs **must not** contain whitespace characters (spaces, tabs, newlines)
- IDs **must not** be empty
- All characters except whitespace and newlines are valid in IDs

.. note:: For strict PHYLIP format compatibility, use IDs of 10 characters or fewer.

Examples
--------
Lower triangular matrix::

    4
    Seq1
    Seq2	1.5
    Seq3	2.0	1.0
    Seq4	3.0	2.5	1.8

Square matrix::

    4
    Seq1	0.0	1.5	2.0	3.0
    Seq2	1.5	0.0	1.0	2.5
    Seq3	2.0	1.0	0.0	1.8
    Seq4	3.0	2.5	1.8	0.0


    https://phylipweb.github.io/phylip/doc/distance.html
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.stats.distance import DistanceMatrix
from skbio.io import create_format, PhylipFormatError

phylip_dm = create_format("phylip_dm")


@phylip_dm.sniffer()
def _phylip_dm_sniffer(fh):
    try:
        header = next(_line_generator(fh))
        n_seqs = _validate_header(header)

        # Collect the first 3 lines for validation
        lines = []
        for line_no in range(3):
            line = next(_line_generator(fh))
            lines.append(line)

        # Try relaxed format first
        try:
            for line_no, line in enumerate(lines):
                _validate_line(line, n_seqs, line_no, strict=False)
            return True, {}
        except PhylipFormatError:
            pass
        # Try strict format
        try:
            for line_no, line in enumerate(lines):
                _validate_line(line, n_seqs, line_no, strict=True)
            return True, {}
        except PhylipFormatError:
            pass

        return False, {}

    except (StopIteration, PhylipFormatError):
        return False, {}


@phylip_dm.reader(DistanceMatrix)
def _phylip_dm_to_distance_matrix(fh, cls=None, strict=False):
    if cls is None:
        cls = DistanceMatrix
    return _phylip_to_dm(cls, fh, strict=strict)


@phylip_dm.writer(DistanceMatrix)
def _distance_matrix_to_phylip(obj, fh):
    _matrix_to_phylip(obj, fh, delimiter="\t")


def _phylip_to_dm(cls, fh, strict=False):
    data = _parse_phylip_dm_raw(fh, strict=strict)
    dists = [x[0] for x in data]
    ids = [x[1] for x in data]
    # If it's in lower triangular form we convert it into condensed form to pass to
    # DistanceMatrix.
    if len(dists[0]) == 0:
        dists = [
            dists[row][col]
            for col in range(len(dists))
            for row in range(col + 1, len(dists))
        ]
    return cls(np.array(dists, dtype=float), ids)


def _matrix_to_phylip(obj, fh, delimiter):
    n_samples = obj.shape[0]
    if n_samples < 2:
        raise PhylipFormatError(
            "DistanceMatrix can only be written in PHYLIP format if there are at least"
            " two samples in the matrix."
        )
    ids = obj.ids
    fh.write(f"{str(n_samples)}\n")

    # if the DistanceMatrix is in redundant form, it will be written in redundant form
    if not obj._flags["CONDENSED"]:
        for id_, vals in zip(ids, obj.data):
            fh.write(id_)
            fh.write(delimiter)
            fh.write(delimiter.join(np.asarray(vals, dtype=str)))
            fh.write("\n")
    # if the DistanceMatrix is in condensed form, it will be written in lower
    # trianglular form
    else:
        for i, id_ in enumerate(ids):
            fh.write(id_)
            if i > 0:
                fh.write(delimiter)
                fh.write(delimiter.join(np.asarray(obj[id_][:i], dtype=str)))
            fh.write("\n")


def _line_generator(fh):
    """Remove linebreak characters and yield lines."""
    for line in fh:
        yield line.rstrip("\n")


def _validate_header(header):
    header_vals = header.split()
    try:
        (n_seqs,) = [int(x) for x in header_vals]
        if n_seqs < 1:
            raise PhylipFormatError("The number of sequences must be positive.")
    except ValueError:
        raise PhylipFormatError(
            "Found non-header line when attempting to read the 1st record "
            f"(header line should have a single integer): {header}"
        )
    return n_seqs


def _validate_line(line, n_seqs, n_dists, strict=False):
    """Check that each line contains the expected number of values.

    Parameters
    ----------
    line
        The line of text being validated.
    n_seqs
        The number of sequences in the matrix.
    n_dists
        The expected number of distances in the line. When a matrix is square, n_dists
        is equal to n_seqs. When a matrix is lower triangle, n_dists is equal to the
        current line number (minus the header), i.e. the 0th non-header line should have
        0 distances, the 1st non-header line should have 1 distance, the 2nd non-header
        line should have 2 distances.

    """
    if not line:
        raise PhylipFormatError("Empty lines are not allowed.")
    if strict:
        id = line[:10].strip()
        dists = line[10:].split()
    else:
        split_line = line.split()
        # IDs are separated from values by whitespace.
        id = split_line[0]
        dists = split_line[1:]
    # This check handles lower triangle matrices. We expects 0 distances on the first
    # non-header line, a single distance on the second non-header line, and two
    # distances on the third non-header line.
    dists_len = len(dists)
    if dists_len != n_dists:
        # If there are more distances than expected for a lower triangle matrix, we
        # expect that it is a square matrix. In this case we check that the number of
        # distances matches the value specified in the header.
        if dists_len != n_seqs:
            raise PhylipFormatError(
                f"The number of distances {dists} is not {n_seqs} as specified in the "
                "header."
            )

    return (dists, id), dists_len


def _parse_phylip_dm_raw(fh, strict=False):
    """Raw parser for PHYLIP formatted distance matrix files."""
    try:
        header = next(_line_generator(fh))
    except StopIteration:
        raise PhylipFormatError("This file is empty.")
    n_seqs = _validate_header(header)

    n_dists = 0
    data = []
    lengths = []

    for line in _line_generator(fh):
        data_, length = _validate_line(line, n_seqs, n_dists, strict=strict)
        data.append(data_)
        lengths.append(length)
        n_dists += 1

    if len(data) != n_seqs:
        raise PhylipFormatError(
            f"The number of sequences is not {n_seqs} as specified in the header."
        )

    # Ensure that no matrix data was accidentally parsed as IDs.
    # The logic here is that either all the distance arrays should be the same length
    # (square format) or all the distance array lengths should be sequentially
    # increasing (lower triangular). If neither is true, then something is wrong.
    is_square = all(L == lengths[0] for L in lengths)
    is_lower_tri = all(lengths[i] == i for i in range(len(lengths)))

    if not (is_square or is_lower_tri):
        raise PhylipFormatError(
            f"Inconsistent distance counts detected: {lengths}. "
            f"This may indicate that sequence IDs contain numeric characters "
            f"being parsed as distances in strict format. "
            f"Expected either all {n_seqs} (square) or 0,1,2,... "
            f"(lower triangular)."
        )

    return data
