"""PHYLIP distance matrix format (:mod:`skbio.io.format.phylip_dm`)
================================================================

.. currentmodule:: skbio.io.format.phylip_dm

The PHYLIP file format stores pairwise distance matrices. This format is commonly
used in phylogenetic analysis and is compatible with tools in the PHYLIP package.

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

Relaxed vs. Strict PHYLIP
^^^^^^^^^^^^^^^^^^^^^^^^^^
scikit-bio supports both **relaxed** and **strict** PHYLIP formats:

**Relaxed PHYLIP** (default):
    - Sequence IDs can have arbitrary length
    - IDs and distance values are separated by whitespace (spaces or tabs)
    - IDs **must not** contain whitespace
    - This is the default format for both reading and writing

**Strict PHYLIP** (optional):
    - Sequence IDs must be exactly 10 characters (padded or truncated)
    - Characters 1-10 are the ID, remaining characters are distance values
    - IDs **may** contain whitespace (e.g., "Sample 01 ")
    - Enable by setting ``strict=True`` when reading

.. note:: scikit-bio writes in relaxed format by default. For strict format
   compatibility with legacy PHYLIP tools, ensure your IDs are 10 characters or
   fewer and do not contain whitespace.

Sequential vs. Interleaved
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The original PHYLIP specification describes both sequential and interleaved formats.
scikit-bio supports **sequential** format only, where each row of the matrix appears
on a single line.

.. note:: Interleaved PHYLIP formats are not currently supported.

Header Section
^^^^^^^^^^^^^^
The header consists of a single line with a single positive integer (``n``) that
specifies the number of sequences in the matrix. This **must** be the first line
in the file. The integer may be preceded by optional whitespace.

.. note:: scikit-bio writes the PHYLIP format header without preceding spaces.
   Empty lines are not allowed between the header and the matrix.

Matrix Section
^^^^^^^^^^^^^^
The matrix section immediately follows the header. It consists of ``n`` lines (rows),
one for each sequence. Each row consists of a sequence identifier (ID) followed by
the distance values for that sequence, separated by whitespace.

**Square matrices**: Each row contains ``n`` distance values (the full row of the
distance matrix, including the diagonal).

**Lower triangular matrices**: Row ``i`` contains ``i`` distance values (only the
values below the diagonal). The first row contains no distance values, just the ID.

Sequence IDs
^^^^^^^^^^^^
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

Reading PHYLIP Files
--------------------
Use the standard scikit-bio I/O interface to read PHYLIP distance matrices:

Reading a file::

    >>> from skbio import DistanceMatrix
    >>> dm = DistanceMatrix.read(  # doctest: +SKIP
    ...     'distance_matrix.phylip', format='phylip_dm')

Reading with strict format parsing::

    >>> dm = DistanceMatrix.read(  # doctest: +SKIP
    ...     'strict_matrix.phylip', format='phylip_dm', strict=True)

Reading from a file handle::

    >>> with open('distance_matrix.phylip', 'r') as f:  # doctest: +SKIP
    ...     dm = DistanceMatrix.read(f, format='phylip_dm')

Writing PHYLIP Files
--------------------
Use the standard scikit-bio I/O interface to write PHYLIP distance matrices. You can
choose between lower triangular (more compact) or square format.

Writing to lower triangular format (default)::

    >>> dm.write('output.phylip', format='phylip_dm')  # doctest: +SKIP
    >>> # or explicitly:
    >>> dm.write(  # doctest: +SKIP
    ...     'output.phylip', format='phylip_dm', lower_tri=True)

Writing to square format::

    >>> dm.write(  # doctest: +SKIP
    ...     'output_square.phylip', format='phylip_dm', lower_tri=False)

Writing to a file handle::

    >>> with open('output.phylip', 'w') as f:  # doctest: +SKIP
    ...     dm.write(f, format='phylip_dm', lower_tri=True)

.. note:: The choice of output format (lower triangular vs. square) is independent
   of how the DistanceMatrix was created. You can write any DistanceMatrix to either
   format.

Format Examples
---------------
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

Strict format with 10-character IDs::

    4
    Sample_001  0.0 1.5 2.0 3.0
    Sample_002  1.5 0.0 1.0 2.5
    Sample_003  2.0 1.0 0.0 1.8
    Sample_004  3.0 2.5 1.8 0.0

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
    >>> dm.write('output.phylip', format='phylip_dm')  # doctest: +SKIP
    >>> # IDs in output will be: 'Sample_01', 'Sample_02', 'Sample_03'

Common Errors
-------------
**"Inconsistent distance counts detected"**: This error typically occurs when:
    - IDs contain whitespace in relaxed format (use ``strict=True`` if needed)
    - The file has irregular formatting or wrong number of values per row

**"The number of distances is not N as specified in the header"**: This occurs when:
    - A row has too many or too few distance values
    - IDs contain whitespace and are being parsed as distance values

**"Empty lines are not allowed"**: PHYLIP format does not allow blank lines between
the header and matrix or within the matrix itself.

References
----------
.. [1] Felsenstein, J. PHYLIP (Phylogeny Inference Package) version 3.6.
   Distributed by the author. Department of Genome Sciences, University of
   Washington, Seattle. 2005.
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
def _distance_matrix_to_phylip(obj, fh, lower_tri=True):
    _matrix_to_phylip(obj, fh, delimiter="\t", lower_tri=lower_tri)


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


def _matrix_to_phylip(obj, fh, delimiter, lower_tri):
    n_samples = obj.shape[0]
    if n_samples < 2:
        raise PhylipFormatError(
            "DistanceMatrix can only be written in PHYLIP format if there are at least"
            " two samples in the matrix."
        )
    ids = obj.ids
    fh.write(f"{str(n_samples)}\n")

    # Default is to write lower triangle
    if lower_tri:
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
    else:
        for id_, vals in zip(ids, obj.data):
            id_w = _remove_whitespace(id_)
            fh.write(id_w)
            fh.write(delimiter)
            fh.write(delimiter.join(np.asarray(obj[id_], dtype=str)))
            fh.write("\n")


def _remove_whitespace(id):
    """Replace whitepace with underscores in ids."""
    return id.replace(" ", "_")


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
                f"The number of distances ({len(dists)}) is not ({n_seqs}) as "
                f"specified in the header. It may be the case that parsing failed due "
                f"to whitespace in the sequence IDs. The first distance value parsed "
                f"is ({dists[0]}), which should be a float. Whitespace in IDs is only "
                f"supported when the 'strict' parameter is set to True."
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
            f"This may indicate that sequence IDs contain whitespace. IDs may only "
            f"contain whitespace if the strict parameter is set to True. "
            f"Expected either all {n_seqs} (square) or 0,1,2,... "
            f"(lower triangular)."
        )

    return data
