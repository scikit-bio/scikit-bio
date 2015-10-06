"""
BLAST+6 format (:mod:`skbio.io.format.blast6`)
==============================================

.. currentmodule:: skbio.io.format.blast6

The BLAST+6 format (``blast+6``) stores strings, integers, and floats relating
to the similarities between multiple DNA sequences. These values are stored in
a simple tabular format with no labels on individual columns. The values are
separated by the tab character.

An example BLAST+6-formatted file comparing two protein sequences, taken
from [1]_::

    moaC<tab>gi|15800534|ref|NP_286546.1|<tab>100.00<tab>161<tab>0<tab>0<tab>1\
    <tab>161<tab>1<tab>161<tab>3e-114<tab>330
    moaC<tab>gi|170768970|ref|ZP_02903423.1|<tab>99.38<tab>161<tab>1<tab>0\
    <tab>1<tab>161<tab>1<tab>161<tab>9e-114<tab>329

Format Support
--------------
**Has Sniffer: No**
**State: Experimental as of 0.4.0-dev.**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`pandas.DataFrame`                                        |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
BLAST+6 format is a tabular text-based format produced by BLAST+ format 6,
containing data related to passed sequences, tab-separated, with no column
labels or headers of any kind. This format includes both BLAST+ output format
6 and BLAST output format 8.

When blasting sequences using Ubuntu's command line, the user may specify what
data is returned.

For example::

    $ blastp -query queries.fasta -subject subjects.fasta -evalue 0.0001\
    -max_target_seqs 2 -outfmt "6 qseqid sseqid bitscore qstart sstart"
    $ query1	subject2	16.9	1	3

In this case, the blastp command returns custom values, which were specified
by the user to be qseqid, sseqid, bitscore, qstart, and sstart. These names
are retained, but as the returned values show, these names do not mark the
columns. The user can also pass ``-outfmt 6`` instead, which will set the
columns to their default values.

For example::

    $ blastp -query queries.fasta -subject subjects.fasta -evalue 0.0001\
    -max_target_seqs 2 -outfmt 6
    $ query1	subject2	100.00	8	0	0	1	8	3	10	9e-05	16.9

If the user does not pass custom values, the blastp command will use its
default column values.

Types of Data
^^^^^^^^^^^^^
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
|ppos       |Percentage of positive-scoring matc\|float|
|           |hes                                 |     |
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
|staxids    |Unique Subject Taxonomy ID(s), sepa\|str  |
|           |rated by a ';' (in numerical order) |     |
+-----------+------------------------------------+-----+
|sscinames  |Unique Subject Scientific Name(s),  |str  |
|           |separated by a ';'                  |     |
+-----------+------------------------------------+-----+
|scomnames  |Unique Subject Common Name(s), sepa\|str  |
|           |rated by a ';'                      |     |
+-----------+------------------------------------+-----+
|sblastnames|unique Subject Blast Name(s), separ\|str  |
|           |ated by a ';' (in alphabetical      |     |
|           |order)                              |     |
+-----------+------------------------------------+-----+
|sskingdoms |unique Subject Super Kingdom(s), se\|str  |
|           |parated by a ';' (in alphabetical   |     |
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

When a blastp data type returns NA values, blastp displays it as
``'N/A'``. When read into a pd.DataFrame, these values are converted
into ``np.nan``.

List of data and descriptions taken from BLAST Command Line Applications User
Manual [2]_.

Format Parameters
-----------------
The following format parameters are available in ``blast+6`` format only:

-``default_columns``: ``False`` by default. If ``True``, will set the column
 labels to the default BLAST+ column labels, which are qseqid, sseqid, pident,
 length, mismatch, gapopen, qstart, qend, sstart, send, evalue, and bitscore.

-``columns``: ``None`` by default. Accepts a list of labels which will be
 assigned to the returned pd.Dataframe.

Examples
--------
An example of reading in ``blast+6`` format with default column labels.

>>> from io import StringIO
>>> import skbio.io as io
>>> import pandas as pd
>>> fs = '\\n'.join([
...     'moaC\\tgi|15800534|ref|NP_286546.1|\\t100.00\\t161\\t0\\t0\\t1\\t161\
\\t1\\t161\\t3e-114\\t330',
...     'moaC\\tgi|170768970|ref|ZP_02903423.1|\\t99.38\\t161\\t1\\t0\\t1\\t\
161\\t1\\t161\\t9e-114\\t329'
... ])
>>> fh = StringIO(fs)

Load BLAST+6 data into a pd.Dataframe and specify that default columns should
be used.

>>> df = io.read(fh, format="blast+6", into=pd.DataFrame, default_columns=True)
>>> df # doctest: +NORMALIZE_WHITESPACE
  qseqid                           sseqid  pident  length  mismatch  gapopen \\
0   moaC     gi|15800534|ref|NP_286546.1|  100.00     161         0        0
1   moaC  gi|170768970|ref|ZP_02903423.1|   99.38     161         1        0
<BLANKLINE>
   qstart  qend  sstart  send         evalue  bitscore
0       1   161       1   161  3.000000e-114       330
1       1   161       1   161  9.000000e-114       329

An example of reading in ``blast+6`` format with custom column labels.

>>> from io import StringIO
>>> import skbio.io as io
>>> import pandas as pd
>>> fs = '\\n'.join([
...     'moaC\\t100.00\\t0\\t161\\t0\\t161\\t330\\t1',
...     'moaC\\t99.38\\t1\\t161\\t0\\t161\\t329\\t1'
... ])
>>> fh = StringIO(fs)

Load BLAST+6 data into a pd.Dataframe and specify your column types.

>>> df = io.read(fh, format="blast+6", into=pd.DataFrame,
...           columns=['qseqid', 'pident', 'mismatch', 'length', 'gapopen',
...                    'qend', 'bitscore', 'sstart'])
>>> df # doctest: +NORMALIZE_WHITESPACE
  qseqid  pident  mismatch  length  gapopen  qend  bitscore  sstart
0   moaC  100.00         0     161        0   161       330       1
1   moaC   99.38         1     161        0   161       329       1

References
----------
.. [1] http://blastedbio.blogspot.com/2014/11/column-headers-in-blast-tabular-\
and-csv.html
.. [2] http://www.ncbi.nlm.nih.gov/books/NBK279675/
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import pandas as pd

from skbio.io import create_format, BLAST6FormatError

blast6 = create_format('blast+6')

_possible_columns = {'qseqid': str, 'qgi': int, 'qacc': str, 'qaccver': str,
                     'qlen': int, 'sseqid': str, 'sallseqid': str, 'sgi': int,
                     'sallgi': int, 'sacc': str, 'saccver': str,
                     'sallacc': str, 'slen': int, 'qstart': int, 'qend': int,
                     'sstart': int, 'send': int, 'qseq': str, 'sseq': str,
                     'evalue': float, 'bitscore': float, 'score': int,
                     'length': int, 'pident': float, 'nident': int,
                     'mismatch': int, 'positive': int, 'gapopen': int,
                     'gaps': int, 'ppos': float, 'frames': str, 'qframe': int,
                     'sframe': int, 'btop': int, 'staxids': str,
                     'sscinames': str, 'scomnames': str, 'sblastnames': str,
                     'sskingdoms': str, 'stitle': str, 'salltitles': str,
                     'sstrand': str, 'qcovs': int, 'qcovhsp': int}
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
                possible_columns = []
                for key in _possible_columns:
                    possible_columns.append(key)
                raise ValueError("Your column name: %s is not valid."
                                 "The valid column names are:\n%s"
                                 "" % (column, possible_columns))

    lineone = pd.read_csv(fh, na_values='N/A', sep='\t', header=None,
                          keep_default_na=False, nrows=1)

    if len(lineone.columns) != len(columns):
        raise BLAST6FormatError("The specified number of columns: %d does not"
                                " match the number of columns in the file: "
                                "%d." % (len(columns), len(lineone.columns)))

    fh.seek(0)
    df = pd.read_csv(fh, na_values='N/A', sep='\t', header=None,
                     keep_default_na=False, names=columns,
                     dtype=_possible_columns)

    return df
