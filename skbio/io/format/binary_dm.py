"""
Simple binary dissimilarity matrix format (:mod:`skbio.io.format.binary_dm`)
============================================================================

.. currentmodule:: skbio.io.format.binary_dm

The binary dissimlarity matrix format stores numeric square matrix data. The
values contained can be interpreted as dissimilarities or distances between
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
`identifiers` dataset is of a variable length unicode type, while the
`dissimilarities` dataset are floating point. The shape of the `identifiers` is
`(N,)`, and the shape of the `dissimilarities` is `(N, N)`.

The dissimilarity between `identifiers[i]` and `identifiers[j]` is interpreted
to be the value at `dissimiarlities[i, j]`. `i` and `j` are integer indices.

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

from skbio.stats.distance import DissimilarityMatrix, DistanceMatrix


_format_header = {
    'format': "Binary DisSimilarity Matrix",
    'version': '2019.7',
    'url': 'https://github.com/biocore/scikit-bio'
    }

def _h5py_mat_to_skbio_mat(cls, fh):
    pass

def _skbio_mat_to_h5py_mat(cls, fh):
    pass

def _get_header(fh):
    pass

def _parse_ids(fh):
    pass

def _verify_dimensions(fh):
    pass

def _bytes_decoder(x):
    return x.decode('ascii')

def _passthrough_decoder(x):
    return x

def _set_header(h5grp):
    """Set format spec header information"""
    for k, v in _format_header.items():
        h5grp.attrs[k] = v

