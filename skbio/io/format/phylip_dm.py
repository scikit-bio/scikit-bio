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

dtype : str or dtype, optional
    The data type of the underlying matrix data. Only relevant when reading from a
    file. Default is "float64", which maps to ``np.float64``. The only other available
    option is "float32" (or ``np.float32``).

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
    n_objs = None
    square = None
    strict = False

    # First line defines number of objects.
    try:
        n_objs = _validate_header(next(fh))
    except (StopIteration, PhylipDMFormatError):
        return False, {}

    # Remaining lines define matrix body (ID + distance values). We will collect up to
    # 5 lines for validation.
    lines = []
    for _ in range(5):
        try:
            lines.append(next(fh).rstrip())
        except StopIteration:
            break

    # must have at least one object
    if (n_lines := len(lines)) < 1:
        return False, {}

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
            id_, *vals = line.rstrip().split()

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

    # return True, {"strict": strict}
    return True, {"strict": strict, "layout": "square" if square else "lower"}


@phylip_dm.reader(DistanceMatrix)
def _phylip_dm_to_distance_matrix(
    fh, cls=None, strict=False, layout="lower", dtype="float64"
):
    if cls is None:
        cls = DistanceMatrix
    dtype = _check_dtype(dtype)
    square = _check_layout(layout)
    data, ids = _phylip_dm_to_matrix(fh, strict, square, dtype)
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


def _phylip_dm_to_matrix(fh, strict, square, dtype):
    """Parse a PHYLIP formatted distance matrix file.

    Parameters
    ----------
    fh : iterator of str
        File handle of a PHYLIP formatted distance matrix.
    strict : bool
        Strict (True) or relaxed (False) object ID format.
    square : bool
        Square (True) or lower triangular (False) layout of matrix body.
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

    """
    # read header and determine number of objects
    try:
        header = next(fh)
    except StopIteration:
        raise PhylipDMFormatError("This file is empty.")
    n_objs = _validate_header(header)

    # allocate array space
    data = np.empty((n_objs, n_objs), dtype=dtype)
    ids = []

    msg_nrow = "The number of rows is not {} as specified in the header."
    msg_ncol = (
        "The number of columns in line {} is not {} as expected."
        "Expected either all {} (square) or 0,1,2,... (lower triangular)."
    )
    msg_plus = (
        "This may indicate that object IDs contain whitespace. IDs may contain"
        "whitespace only in the strict format."
    )

    for i, line in enumerate(fh):
        line = line.rstrip()
        if not line:
            raise PhylipDMFormatError("Empty lines are not allowed.")

        # IDs are strictly 10 characters or less
        if strict:
            id_ = line[:10].strip()
            rest = line[10:]

        # IDs are separated from values by whitespace.
        else:
            try:
                id_, rest = line.split(None, 1)
            except ValueError:
                id_, rest = line.strip(), ""

        vals = np.fromstring(rest, sep=" ", dtype=dtype)
        n_vals = vals.size

        if square:
            if n_vals != n_objs:
                raise PhylipDMFormatError()
        else:
            if n_vals != i:
                raise PhylipDMFormatError()

        try:
            data[i, :n_vals] = vals
        except IndexError:
            raise PhylipDMFormatError(
                f"The number of objects is not {n_objs} as specified in the header."
            )

        ids.append(id_)

    if len(ids) < n_objs:
        raise PhylipDMFormatError(
            f"The number of objects is not {n_objs} as specified in the header."
        )

    if not square:
        np.fill_diagonal(data, 0.0)
        upper = np.triu_indices(n_objs, k=1)
        data[upper] = data.T[upper]

    return data, ids

    # n_dists = 0
    # lengths = []
    # ids = []
    # for line in fh:
    #     id_, data_, length_ = _validate_line(
    #         line.rstrip(), n_objs, n_dists, strict=strict, dtype=dtype
    #     )
    #     try:
    #         data[n_dists, :length_] = data_
    #     except IndexError:
    #         raise PhylipDMFormatError(
    #             f"The number of objects is not {n_objs} as specified in the header."
    #         )
    #     lengths.append(length_)
    #     ids.append(id_)
    #     n_dists += 1

    # if n_dists != n_objs:
    #     raise PhylipDMFormatError(
    #         f"The number of objects is not {n_objs} as specified in the header."
    #     )

    # # Ensure that no matrix data was accidentally parsed as IDs.
    # # The logic here is that either all the distance arrays should be the same length
    # # (square format) or all the distance array lengths should be sequentially
    # # increasing (lower triangular). If neither is true, then something is wrong.
    # is_square = all(L == lengths[0] for L in lengths)
    # is_lower = all(lengths[i] == i for i in range(len(lengths)))

    # if not (is_square or is_lower):
    #     raise PhylipDMFormatError(
    #         f"Inconsistent distance counts detected: {lengths}. "
    #         "This may indicate that object IDs contain whitespace. IDs may only "
    #         "contain whitespace if the strict parameter is set to True. "
    #         f"Expected either all {n_objs} (square) or 0,1,2,... (lower triangular)."
    #     )

    # if is_lower:
    #     np.fill_diagonal(data, 0.0)
    #     upper = np.triu_indices(n_objs, k=1)
    #     data[upper] = data.T[upper]

    # return data, ids


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


def _validate_line(line, n_objs, n_dists, strict=False, dtype="float64"):
    """Check that each line contains the expected number of values.

    Parameters
    ----------
    line : str
        The line of text being validated.
    n_objs : int
        The number of objects in the matrix.
    n_dists : int
        The expected number of distances in the line. When a matrix is square, n_dists
        is equal to n_objs. When a matrix is lower triangle, n_dists is equal to the
        current line number (minus the header), i.e. the 0th non-header line should have
        0 distances, the 1st non-header line should have 1 distance, the 2nd non-header
        line should have 2 distances.

    Returns
    -------
    str
        Object ID.
    ndarray
        Distance values.
    int
        Number of distance values.

    """
    if not line:
        raise PhylipDMFormatError("Empty lines are not allowed.")

    # IDs are strictly 10 characters or less
    if strict:
        id_ = line[:10].strip()
        dists = line[10:]

    # IDs are separated from values by whitespace.
    else:
        try:
            id_, dists = line.split(None, 1)
        except ValueError:
            id_, dists = line.strip(), ""

    # Convert cell values into numbers.
    # `sep=" "` can handle all contiguous whitespaces, like Python's `str.split()`.
    # if `dists` consists of only contiguous whitespaces, output will be [-1], so it is
    # important to strip it before processing.
    dists = np.fromstring(dists, sep=" ", dtype=dtype)

    # This check handles lower triangle matrices. We expect 0 distance on the first
    # non-header line, a single distance on the second non-header line, and two
    # distances on the third non-header line.
    dists_len = dists.size
    if dists_len != n_dists:
        # If there are more distances than expected for a lower triangle matrix, we
        # expect that it is a square matrix. In this case we check that the number of
        # distances matches the value specified in the header.
        if dists_len != n_objs:
            raise PhylipDMFormatError(
                f"The number of distances ({len(dists)}) is not ({n_objs}) as "
                "specified in the header. It may be the case that parsing failed due "
                "to whitespace in the object IDs. The first distance value parsed is "
                f"({dists[0]}), which should be a float. Whitespace in IDs is only "
                "supported when the 'strict' parameter is set to True."
            )

    return id_, dists, dists_len
