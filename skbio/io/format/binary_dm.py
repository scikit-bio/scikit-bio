"""
Simple binary dissimilarity matrix format (:mod:`skbio.io.format.binary_dm`)
============================================================================

.. currentmodule:: skbio.io.format.binary_dm

The binary dissimlarity matrix format stores numeric square hollow matrix data.
The values contained can be interpreted as dissimilarities or distances between
pairs of samples. The format also stores identifiers for each position along
an axis.

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
HDF5 file. Both datatypes are represented by their own datasets. The
`ids` dataset is of a variable length unicode type, while the
`matrix` dataset are floating point. The shape of the `ids` is
`(N,)`, and the shape of the `dissimilarities` is `(N, N)`. The diagonal of
`matrix` are all zeros.

The dissimilarity between `ids[i]` and `ids[j]` is interpreted
to be the value at `matrix[i, j]`. `i` and `j` are integer indices.

.. note:: This file format is most useful for storing large matrices, or when
   it is desirable to facilitate random access to the matrix data.

Format Parameters
-----------------
The only supported format parameter is ``delimiter``, which defaults to the tab
character (``'\\t'``). ``delimiter`` is used to separate elements in the file
format. ``delimiter`` can be specified as a keyword argument when reading from
or writing to a file.
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import h5py
import numpy as np

from skbio.io import create_format, BinaryFormatSymmetryError
from skbio.stats.distance import DissimilarityMatrix, DistanceMatrix


binary_dm = create_format('binary_dm', encoding='binary')


_format_header = {
    'format': "Binary DisSimilarity Matrix",
    'version': '2019.7',
    'url': 'https://github.com/biocore/scikit-bio'
    }

_vlen_dtype = h5py.special_dtype(vlen=str)


@binary_dm.sniffer()
def _binary_dm_sniffer(fh):
    pass


@binary_dm.reader(DissimilarityMatrix)
def _binary_dm_to_dissimilarity(fh):
    pass


@binary_dm.reader(DistanceMatrix)
def _binary_dm_to_distance(fh):
    pass


@binary_dm.writer(DissimilarityMatrix)
def _dissimilarity_to_binary_dm(obj, fh):
    pass


@binary_dm.writer(DistanceMatrix)
def _distance_to_binary_dm(obj, fh):
    pass


def _h5py_mat_to_skbio_mat(cls, fh):
    if issubclass(cls, DistanceMatrix):
        if not fh.attrs['symmetric']:
            raise BinaryFormatSymmetryError("Matrix is not symmetric, cannot "
                                            "construct a DistanceMatrix")
    # a symmetric DissimilarityMatrix is technically okay
    return cls(fh['matrix'], _parse_ids(fh['ids']))


def _skbio_mat_to_h5py_mat(obj, fh):
    _set_header(fh)

    if isinstance(obj, DistanceMatrix):
        fh.attrs['symmetric'] = True
    else:
        fh.attrs['symmetric'] = False

    ids = fh.create_dataset('ids', shape=(len(obj.ids), ), dtype=_vlen_dtype)
    ids[:] = obj.ids
    fh.create_dataset('matrix', data=obj.data)


def _soft_validate_symmetric(fh):
    test = fh['matrix'][:5, :5]
    return np.allclose(np.tril(test), np.tril(test.T))


def _get_header(fh):
    payload = {k: fh.attrs[k] for k in _format_header if k in fh.attrs}
    if set(payload) != set(_format_header):
        return None
    else:
        return payload


def _parse_ids(ids):
    if isinstance(ids[0], bytes):
        return _bytes_decoder(ids)
    else:
        return _passthrough_decoder(ids)


def _verify_dimensions(fh):
    n = len(fh['ids'])
    return fh['matrix'].shape == (n, n)


def _bytes_decoder(x):
    return [i.decode('utf8') for i in x]


def _passthrough_decoder(x):
    return x


def _set_header(h5grp):
    """Set format spec header information"""
    for k, v in _format_header.items():
        h5grp.attrs[k] = v

