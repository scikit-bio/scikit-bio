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
|Yes   |Yes   |:mod:`skbio.stats.distance.PairwiseMatrix`                     |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.stats.distance.SymmetricMatrix`                    |
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

Required datasets:

+-----------+---------------------+------------------------------+
|Datasets   |Value                |Description                   |
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

Optionally, the file can contain several matrices.
In such a case, the 'matrix' dataset will not exist and we will use
'matrix:0' instead.

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
import numpy as np

from skbio.io import create_format
from skbio.stats.distance import PairwiseMatrix, SymmetricMatrix, DistanceMatrix


binary_dm = create_format("binary_dm", encoding="binary")


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
        mat = f.get("matrix:0")
    if mat is None:
        return False, {}

    n = len(ids)
    if mat.shape != (n, n):
        return False, {}

    return True, {}


@binary_dm.reader(PairwiseMatrix)
def _binary_dm_to_pairwise(fh, cls=None):
    if cls is None:
        cls = PairwiseMatrix
    return _h5py_mat_to_skbio_mat_stream(cls, fh)


@binary_dm.reader(SymmetricMatrix)
def _binary_dm_to_symmetric(fh, cls=None):
    if cls is None:
        cls = SymmetricMatrix
    return _h5py_mat_to_skbio_mat_stream(cls, fh)


@binary_dm.reader(DistanceMatrix)
def _binary_dm_to_distance(fh, cls=None):
    if cls is None:
        cls = DistanceMatrix
    return _h5py_mat_to_skbio_mat_stream(cls, fh)


@binary_dm.writer(PairwiseMatrix)
def _pairwise_to_binary_dm(obj, fh):
    return _skbio_mat_to_h5py_mat_stream(obj, fh)


@binary_dm.writer(SymmetricMatrix)
def _symmetric_to_binary_dm(obj, fh):
    return _skbio_mat_to_h5py_mat_stream(obj, fh)


@binary_dm.writer(DistanceMatrix)
def _distance_to_binary_dm(obj, fh):
    return _skbio_mat_to_h5py_mat_stream(obj, fh)


def _h5py_mat_to_skbio_mat_stream(cls, fh):
    with h5py.File(fh, "r") as f:
        dm = _h5py_mat_to_skbio_mat(cls, f)
    return dm


def _h5py_mat_to_skbio_mat(cls, f):
    mat = f.get("matrix")
    if mat is None:
        mat = f.get("matrix:0")
    mat = np.asarray(mat)
    dm = cls(mat, _parse_ids(f["order"]))
    return dm


def _skbio_mat_to_h5py_mat_stream(obj, fh):
    with h5py.File(fh, "w") as f:
        _skbio_mat_to_h5py_mat(obj, f)


def _skbio_mat_to_h5py_mat(obj, f):
    _set_header(f)

    b_ids = [x.encode("utf-8") for x in obj.ids]
    np_ids = np.array(b_ids)
    f.create_dataset("order", data=np_ids)
    f.create_dataset("matrix", data=obj.data)


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


def _verify_dimensions(f):
    ids = f.get("order")

    mat = f.get("matrix")
    if mat is None:
        mat = f.get("matrix:0")

    if (ids is None) or (mat is None):
        return False
    n = len(ids)
    return mat.shape == (n, n)


def _bytes_decoder(x):
    return [i.decode("utf8") for i in x]


def _passthrough_decoder(x):
    return x


def _set_header(f):
    """Set format spec header information."""
    f.create_dataset("format", data=np.array([b"BDSM"]))
    f.create_dataset("version", data=np.array([b"2020.12"]))
