r"""Labeled square matrix format (:mod:`skbio.io.format.lsmat`)
===========================================================

.. currentmodule:: skbio.io.format.lsmat

The labeled square matrix file format (``lsmat``) stores numeric square
matrix data relating a set of objects along each axis. The format also stores
identifiers (i.e., unique labels) for the objects. The matrix data and
identifiers are stored in delimited text format (e.g., TSV or CSV). This format
supports storing a variety of data types including dissimilarity/distance
matrices, similarity matrices and amino acid substitution matrices.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.stats.distance.DissimilarityMatrix`                |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.stats.distance.DistanceMatrix`                     |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
The labeled square matrix and object identifiers are stored as delimited text.
The first line of the file is the header, which must start with the delimiter,
followed by the IDs for all objects in the matrix. Each of the following lines
must contain an object's ID, followed by a numeric (float or integer) vector
relating the object to all other objects in the matrix. The order of objects is
determined by the IDs in the header.

For example, assume we have a 2x2 distance matrix with IDs ``'a'`` and ``'b'``.
When serialized in this format, the distance matrix might look like::

    <del>a<del>b
    a<del>0.0<del>1.0
    b<del>1.0<del>0.0

where ``<del>`` is the delimiter between elements.

Lines containing only whitespace may occur anywhere throughout the file and are
ignored. Lines starting with ``#`` are treated as comments and are ignored.
Comments may only occur *before* the header.

IDs will have any leading/trailing whitespace removed when they are parsed.

.. note:: This file format is most useful for storing small matrices, or when
   it is desirable to represent the matrix in a human-readable format, or
   easily import the file into another program that supports delimited text
   (e.g., a spreadsheet program). If efficiency is a concern, this format may
   not be the most appropriate choice.

Format Parameters
-----------------
The only supported format parameter is ``delimiter``, which defaults to the tab
character (``'\t'``). ``delimiter`` is used to separate elements in the file
format. Examples include tab (``'\t'``) for TSV format and comma (``','``) for
CSV format. ``delimiter`` can be specified as a keyword argument when reading
from or writing to a file.

A special ``delimiter`` is ``None``, which represents a whitespace of arbitrary
length. This value is useful for reading a fixed-width text file. However, it
cannot be automatically determined, nor can it be specified when writing to a
file.

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import csv

import numpy as np

from skbio.stats.distance import DissimilarityMatrix, DistanceMatrix
from skbio.io import create_format, LSMatFormatError


lsmat = create_format("lsmat")


@lsmat.sniffer()
def _lsmat_sniffer(fh):
    header = _find_header(fh)

    if header is not None:
        try:
            dialect = csv.Sniffer().sniff(header)
            delimiter = dialect.delimiter

            ids = _parse_header(header, delimiter)
            first_id, _ = next(_parse_data(fh, delimiter), (None, None))

            if first_id is not None and first_id == ids[0]:
                return True, {"delimiter": delimiter}
        except (csv.Error, LSMatFormatError):
            pass

    return False, {}


@lsmat.reader(DissimilarityMatrix)
def _lsmat_to_dissimilarity_matrix(fh, delimiter="\t"):
    return _lsmat_to_matrix(DissimilarityMatrix, fh, delimiter)


@lsmat.reader(DistanceMatrix)
def _lsmat_to_distance_matrix(fh, delimiter="\t"):
    return _lsmat_to_matrix(DistanceMatrix, fh, delimiter)


@lsmat.writer(DissimilarityMatrix)
def _dissimilarity_matrix_to_lsmat(obj, fh, delimiter="\t"):
    _matrix_to_lsmat(obj, fh, delimiter)


@lsmat.writer(DistanceMatrix)
def _distance_matrix_to_lsmat(obj, fh, delimiter="\t"):
    _matrix_to_lsmat(obj, fh, delimiter)


def _lsmat_to_matrix(cls, fh, delimiter):
    # We aren't using np.loadtxt because it uses *way* too much memory
    # (e.g, a 2GB matrix eats up 10GB, which then isn't freed after parsing
    # has finished). See:
    # http://mail.scipy.org/pipermail/numpy-tickets/2012-August/006749.html

    # Strategy:
    #   - find the header
    #   - initialize an empty ndarray
    #   - for each row of data in the input file:
    #     - populate the corresponding row in the ndarray with floats

    header = _find_header(fh)
    if header is None:
        raise LSMatFormatError(
            "Could not find a header line containing IDs in the "
            "dissimilarity matrix file. Please verify that the file is "
            "not empty."
        )

    ids = _parse_header(header, delimiter)
    num_ids = len(ids)
    data = np.empty((num_ids, num_ids), dtype=np.float64)

    row_idx = -1
    for row_idx, (row_id, row_data) in enumerate(_parse_data(fh, delimiter)):
        if row_idx >= num_ids:
            # We've hit a nonempty line after we already filled the data
            # matrix. Raise an error because we shouldn't ignore extra data.
            raise LSMatFormatError(
                "Encountered extra row(s) without corresponding IDs in " "the header."
            )

        num_vals = len(row_data)
        if num_vals != num_ids:
            raise LSMatFormatError(
                "There are %d value(s) in row %d, which is not equal to the "
                "number of ID(s) in the header (%d)." % (num_vals, row_idx + 1, num_ids)
            )

        expected_id = ids[row_idx]
        if row_id == expected_id:
            data[row_idx, :] = np.asarray(row_data, dtype=float)
        else:
            raise LSMatFormatError(
                "Encountered mismatched IDs while parsing the "
                "dissimilarity matrix file. Found %r but expected "
                "%r. Please ensure that the IDs match between the "
                "dissimilarity matrix header (first row) and the row "
                "labels (first column)." % (str(row_id), str(expected_id))
            )

    if row_idx != num_ids - 1:
        raise LSMatFormatError(
            "Expected %d row(s) of data, but found %d." % (num_ids, row_idx + 1)
        )

    return cls(data, ids)


def _find_header(fh):
    header = None

    for line in fh:
        stripped_line = line.strip()

        if stripped_line and not stripped_line.startswith("#"):
            # Don't strip the header because the first delimiter might be
            # whitespace (e.g., tab).
            header = line
            break

    return header


def _parse_header(header, delimiter):
    tokens = header.rstrip().split(delimiter)

    if delimiter is not None:
        if tokens[0]:
            raise LSMatFormatError(
                "Header must start with delimiter %r." % str(delimiter)
            )
        tokens = tokens[1:]

    return [e.strip() for e in tokens]


def _parse_data(fh, delimiter):
    for line in fh:
        stripped_line = line.strip()

        if not stripped_line:
            continue

        tokens = line.rstrip().split(delimiter)
        id_ = tokens[0].strip()

        yield id_, tokens[1:]


def _matrix_to_lsmat(obj, fh, delimiter):
    delimiter = "%s" % delimiter
    ids = obj.ids
    fh.write(_format_ids(ids, delimiter))
    fh.write("\n")

    for id_, vals in zip(ids, obj.data):
        fh.write("%s" % id_)
        fh.write(delimiter)
        fh.write(delimiter.join(np.asarray(vals, dtype=str)))
        fh.write("\n")


def _format_ids(ids, delimiter):
    return delimiter.join([""] + list(ids))
