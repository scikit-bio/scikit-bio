"""
BLAST+6 format (:mod:`skbio.io.format.blast6`)
=========================================================================

.. currentmodule:: skbio.io.format.blast6

The BLAST+6 format (``blast6``) stores strings, integers, and floats relating
to the similarities between multiple DNA sequences. These values are stored in
a simple tabular format with no labels on individual columns. The values are
separated by the tab character.

An example BLAST+6-formatted file comparing two protein sequences, taken
from [1]_::

    moaC    gi|15800534|ref|NP_286546.1|    100.00  161 0   0   1   161 1
    161 3e-114  330
    moaC    gi|170768970|ref|ZP_02903423.1| 99.38   161 1   0   1   161 1
    161 9e-114  329

Format Support
--------------
**Has Sniffer: No**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`pandas.DataFrame`                                        |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
BLAST+6 format is a tabular text-based format produced by BLAST+ format 6,
containing data related to passed sequences, tab-separated, with no column
labels or headers of any kind.

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

If the user does not pass custom values, the blastp command will use the
default columns, which are qseqid, sseqid, pident, length, mismatch, gapopen,
qstart, qend, sstart, send, evalue, and bitscore.

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
|sallseqid  |All subject Seq-id(s), separated by\|str  |
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
|frames     |Query and subject frames separated \|str  |
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
|sscinames  |Unique Subject Scientific Name(s), \|str  |
|           |separated by a ';'                  |     |
+-----------+------------------------------------+-----+
|scomnames  |Unique Subject Common Name(s), sepa\|str  |
|           |rated by a ';'                      |     |
+-----------+------------------------------------+-----+
|sblastnames|unique Subject Blast Name(s), separ\|str  |
|           |ated by a ';' (in alphabetical\     |     |
|           |order)                              |     |
+-----------+------------------------------------+-----+
|sskingdoms |unique Subject Super Kingdom(s), se\|str  |
|           |parated by a ';' (in alphabetical\  |     |
|           |order)                              |     |
+-----------+------------------------------------+-----+
|stitle     |Subject Title                       |str  |
+-----------+------------------------------------+-----+
|sstrand    |Subject Strand                      |str  |
+-----------+------------------------------------+-----+
|salltitles |All Subject Title(s), separated by \|str  |
|           |a '<>'                              |     |
+-----------+------------------------------------+-----+
|qcovs      |Query Coverage Per Subject          |int  |
+-----------+------------------------------------+-----+
|qcovhsp    |Query Coverage Per HSP              |int  |
+-----------+------------------------------------+-----+

List of data and descriptions taken from BLAST Command Line Applications User
Manual [2]_.

Examples
--------
An example of using _blast6_to_data_frame with default column labels.

>>> from io import StringIO
>>> from skbio import read
>>> import pandas as pd
>>> fs = '\\n'.join([
...     'query1\\tsubject1\\t30.0\\t3\\t5\\t1\\t1\\t3\\t8\\t7\\t0.030\\t13.2',
...     'query2\\tsubject2\\t40.0\\t5\\t4\\t2\\t1\\t6\\t9\\t0\\t0.045\\t15.0'
... ])
>>> fh = StringIO(fs)

Load BLAST+6 data into a DataFrame and specify that default columns should be
used.

>>> df = read(fh, format="blast+6", into=pd.DataFrame, default_columns=True)
>>> df # doctest: +NORMALIZE_WHITESPACE
   qseqid    sseqid  pident  length  mismatch  gapopen  qstart  qend  sstart \\
0  query1  subject1      30       3         5        1       1     3       8
1  query2  subject2      40       5         4        2       1     6       9
<BLANKLINE>
   send  evalue  bitscore
0     7   0.030      13.2
1     0   0.045      15.0

An example of using _blast6_to_data_frame with custom column labels.

>>> from io import StringIO
>>> from skbio import read
>>> import pandas as pd
>>> fs = '\\n'.join([
...     'subject2\\tquery1\\t100\\tN/A\\t0\\tsubject2\\t32\\t100.0',
...     'subject1\\tquery2\\t70\\tN/A\\t0\\tsubject1\\t19\\t70.0',
...     'subject2\\tquery1\\t100\\tN/A\\t0\\tsubject2\\t18\\t50.0'
... ])
>>> fh = StringIO(fs)

Load BLAST+6 data into a DataFrame and specify your column types.

>>> df = read(fh, format="blast+6", into=pd.DataFrame,
...           columns=['sacc', 'qaccver', 'qcovs', 'sblastnames', 'gapopen',
...                    'sallacc', 'score', 'ppos'])
>>> df # doctest: +NORMALIZE_WHITESPACE
       sacc qaccver  qcovs  sblastnames  gapopen   sallacc  score  ppos
0  subject2  query1    100          NaN        0  subject2     32   100
1  subject1  query2     70          NaN        0  subject1     19    70
2  subject2  query1    100          NaN        0  subject2     18    50

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
from skbio.io import create_format
import pandas as pd
import numpy as np

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


@blast6.reader(pd.DataFrame)
def _blast6_to_data_frame(fh, columns=None, default_columns=False):
    if default_columns is True and columns is not None:
        raise ValueError("You cannot use the default columns and also create"
                         " your own columns")
    if default_columns is False and columns is None:
        raise ValueError("You must specify whether default columns or custom"
                         " columns are being used.")
    if columns is not None:
        names = columns
    else:
        names = _default_columns

    df = pd.read_csv(fh, na_values='N/A', sep='\t', header=None)

    if len(list(df.columns.values)) != len(names):
        raise ValueError("The amount of data you have supplied does not match"
                         " the number of columns.")

    try:
        for i in range(len(names)):
            df[i] = df[i].astype(_possible_columns[names[i]])
            df.replace('nan', np.nan, inplace=True)
        df.columns = names
    except ValueError as e:
        raise ValueError('Values have been assigned to incorrect columns. '
                         'Original Pandas error message: '
                         '"%s"' % e)
    return df
