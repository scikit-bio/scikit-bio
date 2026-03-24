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
|Yes   |Yes   |:mod:`skbio.stats.distance.PairwiseMatrix`                     |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.stats.distance.SymmetricMatrix`                    |
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
delimiter : str, optional
    A string used to separate elements in the file format. Can be specified when
    reading from or writing to a file. Default is the tab character (``'\t'``), and
    the format is usually referred to as tab-separated values (TSV). Another common
    choice is comma (``','``), used in the comma-separated values (CSV) format. A
    special ``delimiter`` is ``None``, which represents a whitespace of arbitrary
    length. This value is useful for reading a fixed-width text file. However, it
    cannot be automatically determined, nor can it be specified when writing to a
    file.

dtype : str or dtype, optional
    The data type of the underlying matrix data. Only relevant when reading from a
    file. Default is "float64", which maps to ``np.float64``. The only other available
    option is "float32" (or ``np.float32``).

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

from skbio.stats.distance import PairwiseMatrix, DistanceMatrix, SymmetricMatrix
from skbio.io import create_format, LSMatFormatError


lsmat = create_format("lsmat")


@lsmat.sniffer()
def _lsmat_sniffer(fh):
    header = _find_header(fh)

    if header is not None:
        try:
            dialect = csv.Sniffer().sniff(header)
            delimiter = dialect.delimiter

            # csv delimiter " " is equivalent to Python's str.split(sep=None)
            if delimiter == " ":
                delimiter = None

            ids = _parse_header(header, delimiter)
            first_id, _ = next(_parse_data(fh, delimiter), (None, None))

            if first_id is not None and first_id == ids[0]:
                return True, {"delimiter": delimiter}
        except (csv.Error, LSMatFormatError):
            pass

    return False, {}


@lsmat.reader(PairwiseMatrix)
def _lsmat_to_pairwise_matrix(fh, cls=None, delimiter="\t", dtype="float64"):
    if cls is None:
        cls = PairwiseMatrix
    return _lsmat_to_matrix(cls, fh, delimiter, dtype)


@lsmat.reader(SymmetricMatrix)
def _lsmat_to_symmetric_matrix(fh, cls=None, delimiter="\t", dtype="float64"):
    if cls is None:
        cls = SymmetricMatrix
    return _lsmat_to_matrix(cls, fh, delimiter, dtype)


@lsmat.reader(DistanceMatrix)
def _lsmat_to_distance_matrix(fh, cls=None, delimiter="\t", dtype="float64"):
    if cls is None:
        cls = DistanceMatrix
    return _lsmat_to_matrix(cls, fh, delimiter, dtype)


@lsmat.writer(PairwiseMatrix)
def _pairwise_matrix_to_lsmat(obj, fh, delimiter="\t"):
    _matrix_to_lsmat(obj, fh, delimiter)


@lsmat.writer(SymmetricMatrix)
def _symmetric_matrix_to_lsmat(obj, fh, delimiter="\t"):
    _matrix_to_lsmat(obj, fh, delimiter)


@lsmat.writer(DistanceMatrix)
def _distance_matrix_to_lsmat(obj, fh, delimiter="\t"):
    _matrix_to_lsmat(obj, fh, delimiter)


def _lsmat_to_matrix(cls, fh, delimiter, dtype):
    # We aren't using np.loadtxt because it uses *way* too much memory
    # (e.g, a 2GB matrix eats up 10GB, which then isn't freed after parsing
    # has finished). See:
    # http://mail.scipy.org/pipermail/numpy-tickets/2012-August/006749.html

    # Strategy:
    #   - find the header
    #   - initialize an empty ndarray
    #   - for each row of data in the input file:
    #     - populate the corresponding row in the ndarray with floats

    dtype = np.dtype(dtype)
    if dtype not in (np.float64, np.float32):
        raise TypeError(f"{dtype} is not a supported data type.")

    header = _find_header(fh)
    if header is None:
        raise LSMatFormatError(
            "Could not find a header line containing IDs in the matrix file. Please "
            "verify that the file is not empty."
        )

    ids = _parse_header(header, delimiter)
    num_ids = len(ids)
    data = np.empty((num_ids, num_ids), dtype=dtype)

    # The default separator of str.split (consecutive whitespaces) is equivalent to " "
    # in np.fromstring.
    sep = delimiter if delimiter else " "

    row_idx = -1
    for line in fh:
        line = line.strip()
        if not line:
            continue

        row_idx += 1
        if row_idx >= num_ids:
            # We've hit a nonempty line after we already filled the data matrix. Raise
            # an error because we shouldn't ignore extra data.
            raise LSMatFormatError(
                "Encountered extra row(s) without corresponding IDs in the header."
            )

        row_id, row_data = line.split(delimiter, 1)
        if row_id.rstrip() != ids[row_idx]:
            raise LSMatFormatError(
                "Encountered mismatched IDs while parsing the matrix file. Found %r "
                "but expected %r. Please ensure that the IDs match between the matrix "
                "header (first row) and the row labels (first column)."
                % (row_id.rstrip(), ids[row_idx])
            )

        # np.fromstring is more efficient than str.split then float
        row_data = np.fromstring(row_data, dtype=dtype, sep=sep)

        # The code `data[row_idx, :] = row_data` will raise a ValueError if length
        # doesn't match. However there is an exception: when row_data contains only one
        # element, it will be broadcasted to the entire data[row_idx, :]. Therefore, an
        # explicit check appears to be necessary.
        if row_data.size != num_ids:
            raise LSMatFormatError(
                "There are %d value(s) in row %d, which is not equal to the number of "
                "ID(s) in the header (%d)." % (row_data.size, row_idx + 1, num_ids)
            )
        data[row_idx, :] = row_data

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
        if line := line.rstrip():
            tokens = line.split(delimiter)
            yield tokens[0].strip(), tokens[1:]


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
