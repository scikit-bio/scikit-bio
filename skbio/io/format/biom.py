r"""BIOM-Format (:mod:`skbio.io.format.biom`)
============================================

.. currentmodule:: skbio.io.format.biom

The BIOM-Format (format v2.1.0) is an HDF5-based format to represent sample/feature
counts or relative abundances. It is designed specifically for sparse data.
Internally, it stores the data in both compressed sparse row, and compressed
sparse column representation. It additionally has support for representing sample
and feature metadata.

.. note::

   Internally, BIOM describes features as observations, whereas scikit-bio uses the
   term features.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+-------------------------------------------------------+
|Reader|Writer|                      Object Class                     |
+======+======+=======================================================+
|Yes   |Yes   |:mod:`skbio.table.Table`                               |
+------+------+-------------------------------------------------------+

Format Specification
--------------------
The official format specification for BIOM-Format can be found at [1]_.

Examples
--------
Here we will write an existing BIOM table, and re-read it. Note that the Table
from ``biom`` implicitly gets the ``.write`` method from the IO registry. This
``ByteIO`` object can be a file path in a regular use case.

>>> import io, skbio
>>> f = io.BytesIO()
>>> skbio.table.example_table.write(f)  # doctest: +ELLIPSIS
<_io.BytesIO object at ...>
>>> roundtrip = skbio.read(f, into=skbio.Table)
>>> roundtrip
2 x 3 <class 'biom.table.Table'> with 5 nonzero entries (83% dense)

References
----------
.. [1] http://biom-format.org/documentation/format_versions/biom-2.1.html

"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import h5py

import skbio
from skbio.io import create_format
from skbio.table import Table

from .. import BIOMFormatError


biom = create_format("biom", encoding="binary")


@biom.sniffer()
def _biom_sniffer(fh):
    # this can be buffered, in which case .peek will return the buffer
    # so slice just in case
    magic = fh.peek(8)[:8]

    # From https://en.wikipedia.org/wiki/Hierarchical_Data_Format
    # Note that Wikipedia specifies: "\211HDF\r\n\032\n" which is an ordinal form:
    # >>> ord('\211')
    # 137
    # >>> ord('\x89')
    # 137
    # >>> ord('\032')
    # 26
    # >>> ord('\x1a')
    # 26
    if magic == b"\x89HDF\r\n\x1a\n":
        fp = h5py.File(fh, "r")
        url = fp.attrs.get("format-url")
        version = fp.attrs.get("format-version")

        if url is None or version is None:
            return False, {}
        if url != "http://biom-format.org":
            return False, {}
        if list(version) != [2, 1]:
            return False, {}

        return True, {}
    else:
        return False, {}


@biom.reader(Table)
def _biom_to_table_into(fh):
    return _biom_to_table(fh)


@biom.reader(None)
def _biom_to_table_default(fh):
    # skbio.read('foo.biom', format='biom')
    # will return a generator, that subsequently iterates the table.
    # returning a single item tuple yields expected behavior such that:
    # next(skbio.read('foo.biom', format='biom')) == Table
    return (_biom_to_table(fh),)


def _biom_to_table(fh):
    h5grp = h5py.File(fh, "r")
    return Table.from_hdf5(h5grp)


@biom.writer(Table)
def _sktable_to_biom(obj, fh):
    _table_to_biom(obj, fh)


def _table_to_biom(obj, fh):
    h5grp = h5py.File(fh, "w")
    obj.to_hdf5(h5grp, f"Written by scikit-bio version {skbio.__version__}")
