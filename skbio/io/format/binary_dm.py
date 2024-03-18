"""Simple binary dissimilarity matrix format (:mod:`skbio.io.format.binary_dm`)
============================================================================

.. currentmodule:: skbio.io.format.binary_dm

The Binary DisSimilarity Matrix format (``binary_dm``) encodes a binary
representation for dissimilarity and distance matrices. The format is
designed to facilitate rapid random access to individual rows or columns of
a hollow matrix.

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
The binary dissimilarity matrix and object identifiers are stored within an
HDF5 [1]_ file. Both datatypes are represented by their own datasets. The
`ids` dataset is of a variable length unicode type, while the
`matrix` dataset are floating point. The shape of the `ids` is
`(N,)`, and the shape of the `dissimilarities` is `(N, N)`. The diagonal of
`matrix` are all zeros.

The dissimilarity between `ids[i]` and `ids[j]` is interpreted
to be the value at `matrix[i, j]`. `i` and `j` are integer indices.

Required attributes:

+-----------+---------------------+------------------------------+
|Attribute  |Value                |Description                   |
|           |type                 |                              |
+===========+=====================+==============================+
|format     |string               |A string identifying the file |
|           |                     |as Binary DM format           |
+-----------+---------------------+------------------------------+
|version    |string               |The version of the current    |
|           |                     |Binary DM format              |
+-----------+---------------------+------------------------------+
|matrix     |float32 or float64   |A (N, N) dataset containing   |
|           |                     |the values of the             |
|           |                     |dissimilarity matrix          |
+-----------+---------------------+------------------------------+
|order      |string               |A (N,) dataset of the sample  |
|           |                     |IDs, where N is the total     |
|           |                     |number of IDs                 |
+-----------+---------------------+------------------------------+

.. note:: This file format is most useful for storing large matrices that do
   not need to be represented in a human-readable format. This format is
   especially appropriate for facilitating random access to entries in the
   distance matrix, such as when calculating within and between distances for a
   subset of samples in a large matrix.

References
----------
.. [1] http://www.hdfgroup.org/


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import h5py

from skbio.io import create_format
from skbio.stats.distance import DissimilarityMatrix, DistanceMatrix


binary_dm = create_format("binary_dm", encoding="binary")
_vlen_dtype = h5py.special_dtype(vlen=str)


@binary_dm.sniffer()
def _binary_dm_sniffer(fh):
    try:
        f = h5py.File(fh, "r")
    except OSError:
        return False, {}

    header = _get_header(f)
    if header is None:
        return False, {}

    ids = f.get("order")
    if ids is None:
        return False, {}

    mat = f.get("matrix")
    if mat is None:
        return False, {}

    n = len(ids)
    if mat.shape != (n, n):
        return False, {}

    return True, {}


@binary_dm.reader(DissimilarityMatrix)
def _binary_dm_to_dissimilarity(fh):
    return _h5py_mat_to_skbio_mat(fh)


@binary_dm.reader(DistanceMatrix)
def _binary_dm_to_distance(fh):
    return _h5py_mat_to_skbio_mat(fh)


@binary_dm.writer(DissimilarityMatrix)
def _dissimilarity_to_binary_dm(obj, fh):
    return _skbio_mat_to_h5py_mat(fh)


@binary_dm.writer(DistanceMatrix)
def _distance_to_binary_dm(obj, fh):
    return _skbio_mat_to_h5py_mat(fh)


def _h5py_mat_to_skbio_mat(cls, fh):
    return cls(fh["matrix"], _parse_ids(fh["order"]))


def _skbio_mat_to_h5py_mat(obj, fh):
    _set_header(fh)

    ids = fh.create_dataset("order", shape=(len(obj.ids),), dtype=_vlen_dtype)
    ids[:] = obj.ids
    fh.create_dataset("matrix", data=obj.data)


def _get_header(fh):
    format_ = fh.get("format")
    version = fh.get("version")
    if format is None or version is None:
        return None
    else:
        return {"format": format_[0], "version": version[0]}


def _parse_ids(ids):
    if isinstance(ids[0], bytes):
        return _bytes_decoder(ids)
    else:
        return _passthrough_decoder(ids)


def _verify_dimensions(fh):
    if "order" not in fh or "matrix" not in fh:
        return False
    n = len(fh["order"])
    return fh["matrix"].shape == (n, n)


def _bytes_decoder(x):
    return [i.decode("utf8") for i in x]


def _passthrough_decoder(x):
    return x


def _set_header(h5grp):
    """Set format spec header information."""
    h5grp["format"] = [
        b"BDSM",
    ]
    h5grp["version"] = [
        b"2020.06",
    ]
