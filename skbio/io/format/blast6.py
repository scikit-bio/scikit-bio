"""
BLAST+6 format (:mod:`skbio.io.format.blast6`)
==============================================

.. currentmodule:: skbio.io.format.blast6

The BLAST+6 format (``blast+6``) stores the results of a BLAST [1]_ database
search. The results are stored in a simple tabular format with no column
headers. Values are separated by the tab character.

An example BLAST+6-formatted file comparing two protein sequences, taken
from [2]_ (tab characters represented by ``<tab>``)::

    moaC<tab>gi|15800534|ref|NP_286546.1|<tab>100.00<tab>161<tab>0<tab>0<tab>1\
<tab>161<tab>1<tab>161<tab>3e-114<tab>330
    moaC<tab>gi|170768970|ref|ZP_02903423.1|<tab>99.38<tab>161<tab>1<tab>0\
<tab>1<tab>161<tab>1<tab>161<tab>9e-114<tab>329

Format Support
--------------
**Has Sniffer: No**

**State: Experimental as of 0.4.1.**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`pandas.DataFrame`                                        |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
BLAST+6 format is a tabular text-based format produced by both BLAST+ output
format 6 (``-outfmt 6``) and legacy BLAST output format 8 (``-m 8``). It is
tab-separated and has no column headers. With BLAST+, users can specify the
columns that are present in their BLAST output file by specifying column names
(e.g., ``-outfmt "6 qseqid sseqid bitscore qstart sstart"``), if the default
columns output by BLAST are not desired.

BLAST Column Types
^^^^^^^^^^^^^^^^^^
The following column types are output by BLAST and supported by scikit-bio.
This information is taken from [3]_.

+-----------+------------------------------------+-----+
|Name       |Description                         |Type |
+===========+====================================+=====+
|qseqid     |Query Seq-id                        |str  |
+-----------+------------------------------------+-----+
|qgi        |Query GI                            |int  |
+-----------+------------------------------------+-----+
|qacc       |Query accesion                      |str  |
+-----------+------------------------------------+-----+
|qaccver    |Query accesion.version              |str  |
+-----------+------------------------------------+-----+
|qlen       |Query sequence length               |int  |
+-----------+------------------------------------+-----+
|sseqid     |Subject Seq-id                      |str  |
+-----------+------------------------------------+-----+
|sallseqid  |All subject Seq-id(s), separated by |str  |
|           |a ';'                               |     |
+-----------+------------------------------------+-----+
|sgi        |Subject GI                          |int  |
+-----------+------------------------------------+-----+
|sallgi     |All subject GIs                     |int  |
+-----------+------------------------------------+-----+
|sacc       |Subject accesion                    |str  |
+-----------+------------------------------------+-----+
|saccver    |Subject accesion.version            |str  |
+-----------+------------------------------------+-----+
|sallacc    |All subject accesions               |str  |
+-----------+------------------------------------+-----+
|slen       |Subject sequence length             |int  |
+-----------+------------------------------------+-----+
|qstart     |Start of alignment in query         |int  |
+-----------+------------------------------------+-----+
|qend       |End of alignment in query           |int  |
+-----------+------------------------------------+-----+
|sstart     |Start of alignment in subject       |int  |
+-----------+------------------------------------+-----+
|send       |End of alignment in subject         |int  |
+-----------+------------------------------------+-----+
|qseq       |Aligned part of query sequence      |str  |
+-----------+------------------------------------+-----+
|sseq       |Aligned part of subject sequence    |str  |
+-----------+------------------------------------+-----+
|evalue     |Expect value                        |float|
+-----------+------------------------------------+-----+
|bitscore   |Bit score                           |float|
+-----------+------------------------------------+-----+
|score      |Raw score                           |int  |
+-----------+------------------------------------+-----+
|length     |Alignment length                    |int  |
+-----------+------------------------------------+-----+
|pident     |Percent of identical matches        |float|
+-----------+------------------------------------+-----+
|nident     |Number of identical matches         |int  |
+-----------+------------------------------------+-----+
|mismatch   |Number of mismatches                |int  |
+-----------+------------------------------------+-----+
|positive   |Number of positive-scoring matches  |int  |
+-----------+------------------------------------+-----+
|gapopen    |Number of gap openings              |int  |
+-----------+------------------------------------+-----+
|gaps       |Total number of gaps                |int  |
+-----------+------------------------------------+-----+
|ppos       |Percentage of positive-scoring      |float|
|           |matches                             |     |
+-----------+------------------------------------+-----+
|frames     |Query and subject frames separated  |str  |
|           |by a '/'                            |     |
+-----------+------------------------------------+-----+
|qframe     |Query frame                         |int  |
+-----------+------------------------------------+-----+
|sframe     |Subject frame                       |int  |
+-----------+------------------------------------+-----+
|btop       |Blast traceback operations (BTOP)   |int  |
+-----------+------------------------------------+-----+
|staxids    |Unique Subject Taxonomy ID(s),      |str  |
|           |separated by a ';' (in numerical    |     |
|           |order).                             |     |
+-----------+------------------------------------+-----+
|sscinames  |Unique Subject Scientific Name(s),  |str  |
|           |separated by a ';'                  |     |
+-----------+------------------------------------+-----+
|scomnames  |Unique Subject Common Name(s),      |str  |
|           |separated by a ';'                  |     |
+-----------+------------------------------------+-----+
|sblastnames|unique Subject Blast Name(s),       |str  |
|           |separated by a ';' (in alphabetical |     |
|           |order)                              |     |
+-----------+------------------------------------+-----+
|sskingdoms |unique Subject Super Kingdom(s),    |str  |
|           |separated by a ';' (in alphabetical |     |
|           |order)                              |     |
+-----------+------------------------------------+-----+
|stitle     |Subject Title                       |str  |
+-----------+------------------------------------+-----+
|sstrand    |Subject Strand                      |str  |
+-----------+------------------------------------+-----+
|salltitles |All Subject Title(s), separated by  |str  |
|           |a '<>'                              |     |
+-----------+------------------------------------+-----+
|qcovs      |Query Coverage Per Subject          |int  |
+-----------+------------------------------------+-----+
|qcovhsp    |Query Coverage Per HSP              |int  |
+-----------+------------------------------------+-----+

.. note:: When a BLAST+6-formatted file contains ``N/A`` values, scikit-bio
   will convert these values into ``np.nan``, matching pandas' convention for
   representing missing data.

.. note:: scikit-bio stores columns of type ``int`` as type ``float`` in the
   returned ``pd.DataFrame``. This is necessary in order to allow ``N/A``
   values in integer columns (this is currently a limitation of pandas).

Format Parameters
-----------------
The following format parameters are available in ``blast+6`` format:

- ``default_columns``: ``False`` by default. If ``True``, will use the default
  columns output by BLAST, which are qseqid, sseqid, pident, length, mismatch,
  gapopen, qstart, qend, sstart, send, evalue, and bitscore.

  .. warning::  When reading legacy BLAST files, you must pass
     ``default_columns=True`` because legacy BLAST does not allow users to
     specify which columns are present in the output file.

- ``columns``: ``None`` by default. If provided, must be a list of column names
  in the order they will appear in the file.

.. note:: Either ``default_columns`` or ``columns`` must be provided, as
   ``blast+6`` does not contain column headers.

Examples
--------
Suppose we have a ``blast+6`` file with default columns:

>>> from io import StringIO
>>> import skbio.io
>>> import pandas as pd
>>> fs = '\\n'.join([
...     'moaC\\tgi|15800534|ref|NP_286546.1|\\t100.00\\t161\\t0\\t0\\t1\\t161\
\\t1\\t161\\t3e-114\\t330',
...     'moaC\\tgi|170768970|ref|ZP_02903423.1|\\t99.38\\t161\\t1\\t0\\t1\\t\
161\\t1\\t161\\t9e-114\\t329'
... ])
>>> fh = StringIO(fs)

Read the file into a ``pd.DataFrame`` and specify that default columns should
be used:

>>> df = skbio.io.read(fh, format="blast+6", into=pd.DataFrame,
...                    default_columns=True)
>>> df # doctest: +NORMALIZE_WHITESPACE
  qseqid                           sseqid  pident  length  mismatch  gapopen \\
0   moaC     gi|15800534|ref|NP_286546.1|  100.00   161.0       0.0      0.0
1   moaC  gi|170768970|ref|ZP_02903423.1|   99.38   161.0       1.0      0.0
<BLANKLINE>
   qstart   qend  sstart   send         evalue  bitscore
0     1.0  161.0     1.0  161.0  3.000000e-114     330.0
1     1.0  161.0     1.0  161.0  9.000000e-114     329.0

Suppose we have a ``blast+6`` file with user-supplied (non-default) columns:

>>> from io import StringIO
>>> import skbio.io
>>> import pandas as pd
>>> fs = '\\n'.join([
...     'moaC\\t100.00\\t0\\t161\\t0\\t161\\t330\\t1',
...     'moaC\\t99.38\\t1\\t161\\t0\\t161\\t329\\t1'
... ])
>>> fh = StringIO(fs)

Read the file into a ``pd.DataFrame`` and specify which columns are present
in the file:

>>> df = skbio.io.read(fh, format="blast+6", into=pd.DataFrame,
...                    columns=['qseqid', 'pident', 'mismatch', 'length',
...                             'gapopen', 'qend', 'bitscore', 'sstart'])
>>> df # doctest: +NORMALIZE_WHITESPACE
  qseqid  pident  mismatch  length  gapopen   qend  bitscore  sstart
0   moaC  100.00       0.0   161.0      0.0  161.0     330.0     1.0
1   moaC   99.38       1.0   161.0      0.0  161.0     329.0     1.0

References
----------
.. [1] Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990)
   "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
.. [2] http://blastedbio.blogspot.com/2014/11/column-headers-in-blast-tabular-\
and-csv.html
.. [3] http://www.ncbi.nlm.nih.gov/books/NBK279675/
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from skbio.io import create_format
from skbio.io.format._blast import _parse_blast_data, _possible_columns

blast6 = create_format('blast+6')

_default_columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                    'gapopen', 'qstart', 'qend', 'sstart', 'send',
                    'evalue', 'bitscore']


@blast6.reader(pd.DataFrame, monkey_patch=False)
def _blast6_to_data_frame(fh, columns=None, default_columns=False):
    if default_columns and columns is not None:
        raise ValueError("`columns` and `default_columns` cannot both be"
                         " provided.")
    if not default_columns and columns is None:
        raise ValueError("Either `columns` or `default_columns` must be"
                         " provided.")
    if default_columns:
        columns = _default_columns
    else:
        for column in columns:
            if column not in _possible_columns:
                raise ValueError("Unrecognized column (%r)."
                                 " Supported columns:\n%r" %
                                 (column, set(_possible_columns.keys())))

    return _parse_blast_data(fh, columns, ValueError,
                             "Specified number of columns (%r) does not equal"
                             " number of columns in file (%r).")
