"""

read:

Load dissimilarity matrix from a delimited text file or file path.

Creates a `DissimilarityMatrix` instance from a serialized
dissimilarity matrix stored as delimited text.

`dm_f` can be a file-like or a file path object containing delimited
text. The first line (header) must contain the IDs of each object. The
subsequent lines must contain an ID followed by each dissimilarity
(float) between the current object and all other objects, where the
order of objects is determined by the header line.  For example, a 2x2
dissimilarity matrix with IDs ``'a'`` and ``'b'`` might look like::

    <del>a<del>b
    a<del>0.0<del>1.0
    b<del>1.0<del>0.0

where ``<del>`` is the delimiter between elements.

Parameters
----------
dm_f : iterable of str or str
    Iterable of strings (e.g., open file handle, file-like object, list
    of strings, etc.) or a file path (a string) containing a serialized
    dissimilarity matrix.
delimiter : str, optional
    String delimiting elements in `dm_f`.

Returns
-------
DissimilarityMatrix
    Instance of type `cls` containing the parsed contents of `dm_f`.

Notes
-----
Whitespace-only lines can occur anywhere throughout the "file" and are
ignored. Lines starting with ``#`` are treated as comments and ignored.
These comments can only occur *before* the ID header.

IDs will have any leading/trailing whitespace removed when they are
parsed.

.. note::
    File-like objects passed to this method will not be closed upon the
    completion of the parsing, it is responsibility of the owner of the
    object to perform this operation.

write:

Save the dissimilarity matrix to file in delimited text format.

Parameters
----------
out_f : file-like object or filename
    File-like object to write serialized data to, or name of
    file. If it's a file-like object, it must have a ``write``
    method, and it won't be closed. Else, it is opened and
    closed after writing.
delimiter : str, optional
    Delimiter used to separate elements in output format.

See Also
--------
from_file

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import csv

import numpy as np

from skbio.stats.distance import DissimilarityMatrix, DistanceMatrix
from skbio.io import (register_reader, register_writer, register_sniffer,
                      DMFormatError)

@register_sniffer('dm')
def dm_sniffer(fh):
    valid = False
    kwargs = {}

    header = _find_header(fh)

    if header is None:
        return valid, kwargs

    try:
        dialect = csv.Sniffer().sniff(header)
    except csv.Error:
        return valid, kwargs

    delimiter = dialect.delimiter

    try:
        ids = _parse_header(header, delimiter)
    except DMFormatError:
        return valid, kwargs

    first_id, _ = next(_parse_data(fh, delimiter), (None, None))

    if first_id is not None and first_id == ids[0]:
        valid = True
        kwargs['delimiter'] = delimiter

    return valid, kwargs

@register_reader('dm', DissimilarityMatrix)
def dm_to_DissimilarityMatrix(fh, delimiter='\t'):
    return _dm_to_matrix(DissimilarityMatrix, fh, delimiter)


@register_reader('dm', DistanceMatrix)
def dm_to_DistanceMatrix(fh, delimiter='\t'):
    return _dm_to_matrix(DistanceMatrix, fh, delimiter)


@register_writer('dm', DissimilarityMatrix)
def DissimilarityMatrix_to_dm(obj, fh, delimiter='\t'):
    _matrix_to_dm(obj, fh, delimiter)


@register_writer('dm', DistanceMatrix)
def DistanceMatrix_to_dm(obj, fh, delimiter='\t'):
    _matrix_to_dm(obj, fh, delimiter)


def _dm_to_matrix(cls, fh, delimiter):
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
        raise DMFormatError(
            "Could not find a header line containing IDs in the "
            "dissimilarity matrix file. Please verify that the file is "
            "not empty.")

    ids = _parse_header(header, delimiter)
    num_ids = len(ids)
    data = np.empty((num_ids, num_ids), dtype=np.float64)

    row_idx = -1
    for row_idx, (row_id, row_data) in enumerate(_parse_data(fh, delimiter)):
        if row_idx >= num_ids:
            # We've hit a nonempty line after we already filled the data
            # matrix. Raise an error because we shouldn't ignore extra data.
            raise DMFormatError(
                "Encountered extra row(s) without corresponding IDs in "
                "the header.")

        num_vals = len(row_data)
        if num_vals != num_ids:
            raise DMFormatError(
                "There are %d value(s) in row number %d, which is not "
                "equal to the number of ID(s) in the header (%d)."
                % (num_vals, row_idx + 1, num_ids))

        expected_id = ids[row_idx]
        if row_id == expected_id:
            data[row_idx, :] = np.asarray(row_data, dtype=float)
        else:
            raise DMFormatError(
                "Encountered mismatched IDs while parsing the "
                "dissimilarity matrix file. Found '%s' but expected "
                "'%s'. Please ensure that the IDs match between the "
                "dissimilarity matrix header (first row) and the row "
                "labels (first column)." % (row_id, expected_id))

    if row_idx != num_ids - 1:
        raise DMFormatError(
            "Expected %d row(s) of data, but found %d." % (num_ids, row_idx + 1))

    return cls(data, ids)


def _find_header(fh):
    header = None

    for line in fh:
        stripped_line = line.strip()

        if stripped_line and not stripped_line.startswith('#'):
            # Don't strip the header because the first delimiter might be
            # whitespace (e.g., tab).
            header = line
            break

    return header


def _parse_header(header, delimiter):
    tokens = header.rstrip().split(delimiter)

    if tokens[0]:
        raise DMFormatError(
            "Header must start with delimiter %r." % delimiter)

    return [e.strip() for e in tokens[1:]]


def _parse_data(fh, delimiter):
    for line in fh:
        stripped_line = line.strip()

        if not stripped_line:
            continue

        tokens = line.rstrip().split(delimiter)
        id_ = tokens[0].strip()

        yield id_, tokens[1:]


def _matrix_to_dm(obj, fh, delimiter):
    ids = obj.ids
    fh.write(_format_ids(ids, delimiter))
    fh.write('\n')

    for id_, vals in zip(ids, obj.data):
        fh.write(id_)
        fh.write(delimiter)
        fh.write(delimiter.join(np.asarray(vals, dtype=np.str)))
        fh.write('\n')


def _format_ids(ids, delimiter):
    return delimiter.join([''] + list(ids))
