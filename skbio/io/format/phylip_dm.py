"""PHYLIP distance matrix format (:mod:`skbio.io.format.phylip_dm`)
================================================================

.. currentmodule:: skbio.io.format.phylip_dm

.. versionadded:: 0.7.2

The PHYLIP file format stores pairwise distance matrices. This format is commonly
used in phylogenetic analysis and is compatible with tools in the PHYLIP package.
See [1]_ for the original format description.

An example PHYLIP-formatted distance matrix in lower triangular layout::

    5
    Seq1
    Seq2    1.6866
    Seq3    1.7198  1.5232
    Seq4    1.6606  1.4841  0.7115
    Seq5    1.5243  1.4465  0.5958  0.4631

And its equivalent in square layout::

    5
    Seq1    0.0000  1.6866  1.7198  1.6606  1.5243
    Seq2    1.6866  0.0000  1.5232  1.4841  1.4465
    Seq3    1.7198  1.5232  0.0000  0.7115  0.5958
    Seq4    1.6606  1.4841  0.7115  0.0000  0.4631
    Seq5    1.5243  1.4465  0.5958  0.4631  0.0000


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
describing the number of objects in the matrix, followed by the matrix itself.

Relaxed vs. Strict PHYLIP
^^^^^^^^^^^^^^^^^^^^^^^^^^
scikit-bio supports both **relaxed** and **strict** PHYLIP formats:

**Relaxed PHYLIP** (default):
    - Object IDs can have arbitrary length
    - IDs and distance values are separated by whitespace (spaces or tabs)
    - IDs **must not** contain whitespace
    - This is the default format for both reading and writing

**Strict PHYLIP** (optional):
    - Object IDs must be exactly 10 characters (padded or truncated)
    - Characters 1-10 are the ID, remaining characters are distance values
    - IDs **may** contain whitespace (e.g., "Sample 01 ")
    - Enable by setting ``strict=True`` when reading

.. note:: scikit-bio writes in relaxed format by default. For strict format
   compatibility with legacy PHYLIP tools, ensure your IDs are 10 characters or
   fewer and do not contain whitespace.

Header Section
^^^^^^^^^^^^^^
The header consists of a single line with a single positive integer (``n``) that
specifies the number of objects in the matrix. This **must** be the first line
in the file. The integer may be preceded by optional whitespace.

.. note:: scikit-bio writes the PHYLIP format header without preceding spaces.
   Empty lines are not allowed between the header and the matrix.

Matrix Section
^^^^^^^^^^^^^^
The matrix section immediately follows the header. It consists of ``n`` lines (rows),
one for each object. Each row consists of an object identifier (ID) followed by
the distance values for that object, separated by whitespace. Two alternative
layouts of the matrix body are supported:

**Square matrices**: Each row contains ``n`` distance values (the full row of the
distance matrix, including the diagonal).

**Lower triangular matrices**: Row ``i`` contains ``i`` distance values (only the
values below the diagonal). The first row contains no distance values, just the ID.

.. note:: The original PHYLIP format also defines upper triangular matrices, although
   they are less common and currently not supported by scikit-bio.

Object IDs
^^^^^^^^^^
**Relaxed format** (default):
    - IDs can have arbitrary length
    - IDs **must not** contain whitespace characters (spaces, tabs, newlines)
    - IDs **must not** be empty
    - All characters except whitespace and newlines are valid

**Strict format** (``strict=True``):
    - IDs occupy exactly the first 10 characters of each line
    - IDs **may** contain whitespace
    - IDs are automatically padded or truncated to 10 characters

.. note:: When writing, any whitespace in IDs is automatically replaced with
   underscores to ensure compatibility with relaxed format.


Format Parameters
-----------------

Reader-specific Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^
strict : bool, optional
    Whether the object IDs are in strict (True) or relaxed (False, default) format.

Writer-specific Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^
layout : {'lower', 'square'}, optional
    Layout of the matrix body. Options are "lower" (lower triangle, default) and
    "square" (square).


Examples
--------

Reading PHYLIP Files
^^^^^^^^^^^^^^^^^^^^

Read a PHYLIP distance matrix file into a ``DistanceMatrix`` object:

>>> from skbio import DistanceMatrix
>>> dm = DistanceMatrix.read('input.phy', format='phylip_dm')  # doctest: +SKIP


Read with strict ID format parsing:

>>> dm = DistanceMatrix.read(  # doctest: +SKIP
...     'input_strict.phy', format='phylip_dm', strict=True)

Read from a file handle:

>>> with open('input.phy', 'r') as f:  # doctest: +SKIP
...     dm = DistanceMatrix.read(f, format='phylip_dm')

Writing PHYLIP Files
^^^^^^^^^^^^^^^^^^^^

Use the standard scikit-bio I/O interface to write PHYLIP distance matrices. You can
choose between lower triangular (more compact) or square layout.

Write to lower triangular layout (default):

>>> dm.write('output.phy', format='phylip_dm')  # doctest: +SKIP
>>> # or explicitly:
>>> dm.write('output.phy', format='phylip_dm', layout='lower')  # doctest: +SKIP

Write to square layout:

>>> dm.write('output_square.phy', format='phylip_dm', layout='square')  # doctest: +SKIP

Write to a file handle:

>>> with open('output.phy', 'w') as f:  # doctest: +SKIP
...     dm.write(f, format='phylip_dm')

.. note:: The choice of output format (lower triangular vs. square) is independent
   of how the DistanceMatrix was created. You can write any DistanceMatrix to either
   format.


Working with Whitespace in IDs
-------------------------------
**Relaxed format** does not support whitespace in IDs. If your IDs contain
whitespace and you're using relaxed format, the parser will treat the first
whitespace-delimited token as the ID and subsequent tokens as distance values,
which will cause parsing errors.

**Strict format** supports whitespace in IDs because IDs are positional (first 10
characters). To read files with whitespace in IDs, use ``strict=True``.

When **writing** files, scikit-bio automatically replaces any whitespace in IDs
with underscores to ensure compatibility::

    >>> ids = ['Sample 01', 'Sample 02', 'Sample 03']
    >>> dm = DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]], ids=ids)
    >>> dm.write('output.phy', format='phylip_dm')  # doctest: +SKIP
    >>> # IDs in output will be: 'Sample_01', 'Sample_02', 'Sample_03'

Common Errors
-------------
**"Inconsistent distance counts detected"**: This error typically occurs when:
    - IDs contain whitespace in relaxed format (use ``strict=True`` if needed)
    - The file has irregular formatting or wrong number of values per row

**"The number of distances is not N as specified in the header"**: This occurs when:
    - A row has too many or too few distance values
    - IDs contain whitespace and are being parsed as distance values

**"Empty lines are not allowed"**:
    - PHYLIP format does not allow blank lines between the header and matrix or within
      the matrix itself.

References
----------
.. [1] Felsenstein, J. PHYLIP (Phylogeny Inference Package) version 3.6. Distributed by
   the author. Department of Genome Sciences, University of Washington, Seattle. 2005.
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
from skbio.io import create_format, PhylipDMFormatError

phylip_dm = create_format("phylip_dm")


@phylip_dm.sniffer()
def _phylip_dm_sniffer(fh):
    try:
        header = next(fh).rstrip()
        n_seqs = _validate_header(header)

        # Collect the first 3 lines for validation
        lines = []
        for line_no in range(3):
            line = next(fh).rstrip()
            lines.append(line)

        # Try relaxed format first
        try:
            for line_no, line in enumerate(lines):
                _validate_line(line, n_seqs, line_no, strict=False)
            return True, {}
        except PhylipDMFormatError:
            pass
        # Try strict format
        try:
            for line_no, line in enumerate(lines):
                _validate_line(line, n_seqs, line_no, strict=True)
            return True, {}
        except PhylipDMFormatError:
            pass

        return False, {}

    except (StopIteration, PhylipDMFormatError):
        return False, {}


@phylip_dm.reader(DistanceMatrix)
def _phylip_dm_to_distance_matrix(fh, cls=None, strict=False):
    if cls is None:
        cls = DistanceMatrix
    return _phylip_to_dm(cls, fh, strict=strict)


@phylip_dm.writer(DistanceMatrix)
def _distance_matrix_to_phylip(obj, fh, layout="lower"):
    _matrix_to_phylip(obj, fh, delimiter="\t", layout=layout)


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


def _matrix_to_phylip(obj, fh, delimiter, layout):
    n_samples = obj.shape[0]
    if n_samples < 2:
        raise PhylipDMFormatError(
            "DistanceMatrix can only be written in PHYLIP format if there are at least"
            " two samples in the matrix."
        )
    ids = obj.ids
    fh.write(f"{str(n_samples)}\n")

    # Default is to write lower triangle
    if layout == "lower":
        for i, id_ in enumerate(ids):
            # Here we are replacing whitespace with underscore on write, but we still
            # need to be able to index the DistanceMatrix object by id (which may
            # contain whitespace) so create a separate variable for writing the id
            id_w = _remove_whitespace(id_)
            fh.write(id_w)
            if i > 0:
                fh.write(delimiter)
                fh.write(delimiter.join(np.asarray(obj[id_][:i], dtype=str)))
            fh.write("\n")
    elif layout == "square":
        for id_, vals in zip(ids, obj.data):
            id_w = _remove_whitespace(id_)
            fh.write(id_w)
            fh.write(delimiter)
            fh.write(delimiter.join(np.asarray(obj[id_], dtype=str)))
            fh.write("\n")
    elif layout == "upper":
        raise PhylipDMFormatError("Upper triangular layout is currently not supported.")
    else:
        raise PhylipDMFormatError(f"'{layout}' is not a supported matrix layout.")


def _remove_whitespace(id_):
    """Replace whitepace with underscores in ids."""
    return id_.replace(" ", "_")


def _validate_header(header):
    header_vals = header.split()
    try:
        (n_seqs,) = [int(x) for x in header_vals]
        if n_seqs < 1:
            raise PhylipDMFormatError("The number of objects must be positive.")
    except ValueError:
        raise PhylipDMFormatError(
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
        The number of objects in the matrix.
    n_dists
        The expected number of distances in the line. When a matrix is square, n_dists
        is equal to n_seqs. When a matrix is lower triangle, n_dists is equal to the
        current line number (minus the header), i.e. the 0th non-header line should have
        0 distances, the 1st non-header line should have 1 distance, the 2nd non-header
        line should have 2 distances.

    """
    if not line:
        raise PhylipDMFormatError("Empty lines are not allowed.")
    if strict:
        id_ = line[:10].strip()
        dists = line[10:].split()
    else:
        split_line = line.split()
        # IDs are separated from values by whitespace.
        id_ = split_line[0]
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
            raise PhylipDMFormatError(
                f"The number of distances ({len(dists)}) is not ({n_seqs}) as "
                "specified in the header. It may be the case that parsing failed due "
                "to whitespace in the object IDs. The first distance value parsed is "
                f"({dists[0]}), which should be a float. Whitespace in IDs is only "
                "supported when the 'strict' parameter is set to True."
            )

    return (dists, id_), dists_len


def _parse_phylip_dm_raw(fh, strict=False):
    """Raw parser for PHYLIP formatted distance matrix files."""
    try:
        header = next(fh).rstrip()
    except StopIteration:
        raise PhylipDMFormatError("This file is empty.")
    n_seqs = _validate_header(header)

    n_dists = 0
    data = []
    lengths = []

    for line in fh:
        data_, length = _validate_line(line.rstrip(), n_seqs, n_dists, strict=strict)
        data.append(data_)
        lengths.append(length)
        n_dists += 1

    if len(data) != n_seqs:
        raise PhylipDMFormatError(
            f"The number of objects is not {n_seqs} as specified in the header."
        )

    # Ensure that no matrix data was accidentally parsed as IDs.
    # The logic here is that either all the distance arrays should be the same length
    # (square format) or all the distance array lengths should be sequentially
    # increasing (lower triangular). If neither is true, then something is wrong.
    is_square = all(L == lengths[0] for L in lengths)
    is_lower = all(lengths[i] == i for i in range(len(lengths)))

    if not (is_square or is_lower):
        raise PhylipDMFormatError(
            f"Inconsistent distance counts detected: {lengths}. "
            "This may indicate that object IDs contain whitespace. IDs may only "
            "contain whitespace if the strict parameter is set to True. "
            f"Expected either all {n_seqs} (square) or 0,1,2,... (lower triangular)."
        )

    return data
