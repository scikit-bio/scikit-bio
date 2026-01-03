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

Header Section
^^^^^^^^^^^^^^
The header consists of a single line with a single positive integer (``n``) that
specifies the number of objects in the matrix. This **must** be the first line
in the file. The integer may be preceded by optional whitespace.

.. note:: scikit-bio writes the PHYLIP format header without preceding spaces.

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

.. note:: Empty lines are not allowed between the header and the matrix.

.. note:: The original PHYLIP format also defines upper triangular matrices, although
   they are less common and currently not supported by scikit-bio.

Object IDs
^^^^^^^^^^
scikit-bio supports both **relaxed** and **strict** object ID formats:

**Relaxed format** (default):
    - IDs can have arbitrary length
    - IDs **must not** contain whitespace characters (spaces, tabs, etc.)
    - All characters except whitespace are valid

**Strict format** (``strict=True``):
    - IDs occupy exactly the first 10 characters of each line
    - IDs **may** contain whitespace
    - IDs are automatically padded or truncated to 10 characters

scikit-bio writes in relaxed format by default. For strict format compatibility with
legacy PHYLIP tools, ensure your IDs are 10 characters or fewer and do not contain
whitespace.

.. note:: When writing, any whitespace in IDs is automatically replaced with
   underscores to ensure compatibility with relaxed format.


Format Parameters
-----------------
layout : {'lower', 'square'}, optional
    Layout of the matrix body. Options are "lower" (lower triangle) and "square"
    (square). Applicable to both reading and writing. The layout of the input file
    is automatically inferred during reading, although one can explicitly set this
    parameter to override. Writing defaults to lower triangle.

strict : bool, optional
    Whether the object IDs are in strict (True) or relaxed (False, default) format.
    Only applicable to reading. The ID format of the input file is automatically
    inferred during reading, although one can explicitly set this parameter to
    override. Writing always uses the relaxed format.

dtype : str or dtype, optional
    The data type of the underlying matrix data. Default is "float64", which maps to
    ``np.float64``. The only other available option is "float32" (or ``np.float32``).
    Only relevant when reading from a file.


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
    # First line defines number of objects.
    try:
        n_objs = _parse_header(next(fh))
    except (StopIteration, PhylipDMFormatError):
        return False, {}

    # Remaining lines define matrix body (ID + distance values). We will collect up to
    # 5 lines for validation.
    lines = []
    for _ in range(5):
        try:
            lines.append(next(fh))
        except StopIteration:
            break

    # Must have at least one object.
    if len(lines) < 1:
        return False, {}

    # Sequentially test four combinations of format parameters. If any of these passed,
    # return immediately. If failed, continue to test the next combination.
    # for square in (False, True):
    #     for strict in (False, True):
    #         try:
    #             for i, line in enumerate(lines):
    #                 _parse_line(line, i, n_objs, strict, square, float)
    #         except PhylipDMFormatError:
    #             continue
    #         else:
    #             return True, {
    #                 "layout": "square" if square else "lower",
    #                 "strict": strict,
    #             }

    # return False, {}

    strict = False
    square = None
    n_lines = len(lines)
    i = 0
    while i < n_lines:
        line = lines[i]

        # empty line is prohibited
        if not line:
            return False, {}

        # ID in strict format: up to first 10 characters
        if strict:
            id_ = line[:10].strip()
            vals = line[10:].split()

            # ID cannot be empty
            if not id_:
                return False, {}

        # ID in relaxed format (default): before the first whitespace
        else:
            try:
                id_, *vals = line.rstrip().split()
            except ValueError:
                return False, {}

        n_vals = len(vals)
        to_strict = False

        # Infer layout from the second line: no value => lower triangular layout;
        # same number of values as objects => square layout.
        if i == 0:
            if n_vals == n_objs:
                square = True
            elif n_vals == 0:
                square = False
            elif not strict:
                to_strict = True
            else:
                return False, {}

        # A special case: first ID consumes exactly 10 characters, leaving no
        # whitespace between it and values. As a consequence, there appears to be
        # one less value than objects. When this is observed, we will turn on
        # strict ID format and restart from the second line.

        # For all remaining lines, check if the number of values is expected.
        elif (square and n_vals != n_objs) or (not square and n_vals != i):
            if not strict:
                to_strict = True
            else:
                return False, {}

        # all values must be numeric
        try:
            _ = list(map(float, vals))
        except ValueError:
            if not strict:
                to_strict = True
            else:
                return False, {}

        # Switch to strict ID format and restart the iteration.
        if to_strict:
            strict = True
            square = None
            i = 0
            continue

        i += 1

    return True, {"layout": "square" if square else "lower", "strict": strict}


@phylip_dm.reader(DistanceMatrix)
def _phylip_dm_to_distance_matrix(
    fh, cls=None, layout="lower", strict=False, dtype="float64"
):
    if cls is None:
        cls = DistanceMatrix
    dtype = _check_dtype(dtype)
    square = _check_layout(layout)
    data, ids = _phylip_dm_to_matrix(fh, square, strict, dtype)
    return cls(data, ids)


@phylip_dm.writer(DistanceMatrix)
def _distance_matrix_to_phylip_dm(obj, fh, layout="lower"):
    _matrix_to_phylip_dm(obj, fh, delimiter="\t", layout=layout)


def _check_layout(layout):
    if layout == "lower":
        return False
    elif layout == "square":
        return True
    elif layout == "upper":
        raise PhylipDMFormatError("Upper triangular layout is currently not supported.")
    else:
        raise PhylipDMFormatError(f"'{layout}' is not a supported matrix layout.")


def _check_dtype(dtype):
    if (typ := np.dtype(dtype)) not in (np.float64, np.float32):
        raise TypeError(f"{dtype} is not a supported data type.")
    return typ


def _phylip_dm_to_matrix(fh, square, strict, dtype):
    """Parse a PHYLIP formatted distance matrix file.

    Parameters
    ----------
    fh : iterator of str
        File handle of a PHYLIP formatted distance matrix.
    square : bool
        Square (True) or lower triangular (False) layout of matrix body.
    strict : bool
        Strict (True) or relaxed (False) object ID format.
    dtype : type
        Data type of distance values (float64 or float32).

    Returns
    -------
    data : ndarray of shape (n_objects, n_objects)
        Distance matrix data.
    ids : list of str
        Object IDs.

    Raises
    ------
    PhylipDMFormatError
        If file is empty.
        If matrix body contains more or less rows than defined by the header.

    """
    # parse header and determine number of objects
    try:
        header = next(fh)
    except StopIteration:
        raise PhylipDMFormatError("This file is empty.")
    n_objs = _parse_header(header)

    # allocate array space
    data = np.empty((n_objs, n_objs), dtype=dtype)
    ids = []

    # parse each line of matrix body and append results
    for i, line in enumerate(fh):
        id_, vals, n_vals = _parse_line(line, i, n_objs, square, strict, dtype)
        ids.append(id_)
        try:
            data[i, :n_vals] = vals
        except IndexError:
            raise PhylipDMFormatError(
                f"The number of rows in the matrix body exceeds {n_objs} as specified "
                "in the header."
            )

    if len(ids) < n_objs:
        raise PhylipDMFormatError(
            f"The number of rows in the matrix body ({len(ids)}) does not match the "
            f"number of objects specified in the header ({n_objs})."
        )

    # fill upper triangle of the array
    if not square:
        np.fill_diagonal(data, 0.0)
        upper = np.triu_indices(n_objs, k=1)
        data[upper] = data.T[upper]

    return data, ids


def _matrix_to_phylip_dm(obj, fh, delimiter, layout):
    """Write a distance matrix object into a PHYLIP distance matrix file.

    Parameters
    ----------
    obj : DistanceMatrix
        Distance matrix object.
    fh : file handle
        PHYLIP distance matrix file.
    delimiter : str
        Field delimiter.
    layout : {'square', 'lower'}
        Matrix layout.

    Raises
    ------
    TypeError
        If dtype is invalid.
    PhylipDMFormatError
        If file is empty.

    """
    n_samples = obj.shape[0]
    if n_samples < 2:
        raise PhylipDMFormatError(
            "DistanceMatrix can only be written in PHYLIP format if there are at least"
            " two samples in the matrix."
        )
    ids = obj.ids
    fh.write(f"{str(n_samples)}\n")

    # square layout
    if _check_layout(layout):
        for id_, vals in zip(ids, obj.data):
            id_w = _remove_whitespace(id_)
            fh.write(id_w)
            fh.write(delimiter)
            fh.write(delimiter.join(np.asarray(obj[id_], dtype=str)))
            fh.write("\n")

    # lower triangular layout (default)
    else:
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


def _remove_whitespace(id_):
    """Replace whitepace with underscores in IDs."""
    return id_.replace(" ", "_")


def _parse_header(header):
    """Check that the header line contains a valid number of objects.

    Parameters
    ----------
    header : str
        The header line being validated.

    Returns
    -------
    n_objs
        The number of objects in the matrix.

    Raises
    ------
    PhylipDMFormatError
        If header cannot be converted into an integer.
        If number of objects is zero or negative.

    """
    try:
        n_objs = int(header.strip())
    except ValueError:
        raise PhylipDMFormatError("Header line must be a single integer.")
    if n_objs < 1:
        raise PhylipDMFormatError("The number of objects must be positive.")
    return n_objs


def _parse_line(line, idx, n_objs, square, strict, dtype):
    """Parse each line in the matrix body.

    Parameters
    ----------
    line : str
        Line of text being parsed.
    idx : int
        Line index relative to matrix body (e.g., 2nd line of file has idx = 0).
    n_objs : int
        Number of objects in the matrix.
    square : bool
        Square (True) or lower triangular (False) layout of matrix body.
    strict : bool
        Strict (True) or relaxed (False) object ID format.
    dtype : type
        Data type of distance values (float64 or float32).

    Returns
    -------
    str
        Object ID.
    ndarray
        Distance values.
    int
        Number of distance values.

    Raises
    ------
    PhylipDMFormatError
        If line contains only whitespaces.
        If object ID is empty (strict only).
        If number of objects is zero or negative.

    """
    msg_plus = (
        "This may indicate that object IDs contain whitespace. IDs may contain"
        "whitespace only in the strict format."
    )

    line = line.rstrip()
    if not line:
        raise PhylipDMFormatError("Empty lines are not allowed.")

    # Strict mode: IDs are 10 characters or less.
    if strict:
        id_ = line[:10].strip()
        rest = line[10:]
        if not id_:
            raise PhylipDMFormatError("Empty IDs are not allowed.")

    # Relaxed mode: IDs are separated from values by whitespace.
    else:
        try:
            id_, rest = line.split(None, 1)
        except ValueError:
            id_, rest = line.strip(), ""  # this could happen in the first line

    # Convert cell values into numbers.
    # `sep=" "` can handle all contiguous whitespaces, like Python's `str.split()`.
    # If `rest` consists of only contiguous whitespaces, output will be [-1], so it is
    # important to strip it before processing.
    try:
        vals = np.fromstring(rest, sep=" ", dtype=dtype)
    except ValueError:
        raise PhylipDMFormatError("Non-numeric cell values encountered.")
    else:
        n_vals = vals.size

    # Check if the number of distance values in the line is expected. For square
    # matrices, each line should contain the same number of distances as defined in the
    # header (n_objs).
    if square:
        if n_vals != n_objs:
            raise PhylipDMFormatError(
                f"The number of distance values ({n_vals}) in line {idx + 2} does not "
                f"match the expectation ({n_objs}, as specified in the header)."
            )

    # For lower triangular matrices, we expect 0 distance on the first (non-header)
    # line, a single distance on the second line, and two distances on the third line.
    # So on so forth. Therefore n_vals is equal to the line index (minus the header).
    else:
        if n_vals != idx:
            raise PhylipDMFormatError(
                f"The number of distance values ({n_vals}) in line {idx + 2} does not "
                f"match the expectation ({idx}, which is line number - 2)."
            )

    return id_, vals, n_vals
