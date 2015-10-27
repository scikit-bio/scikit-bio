"""
BLAST+7 format (:mod:`skbio.io.format.blast7`)
==============================================

.. currentmodule:: skbio.io.format.blast7

There are 2 BLAST data format types supported by the reader and sniffer:
BLAST+7 and Legacy BLAST 9. Both are refered to as ``blast+7``. BLAST+7 is
relatively structurally equal to Legacy BLAST format 9, however, this format
has a few differences.

Both formats store the results of a BLAST [1]_ database search. The results
are stored in a simple tabular format with the BLAST version, query, fields,
and database and/or subject specified. Values are separated by the tab
character.

An example BLAST+7-formatted file comparing two protein sequences, taken
from [2]_ (tab characters represented by ``<tab>``):

.. code-block:: none

    # BLASTN 2.2.18+
    # Query: gi|1786181|gb|AE000111.1|AE000111
    # Subject: ecoli
    # Fields: query acc., subject acc., evalue, q. start, q. end, s. st\
art, s. end
    # 5 hits found
    AE000111<tab>AE000111<tab>0.0<tab>1<tab>10596<tab>1<tab>10596
    AE000111<tab>AE000174<tab>8e-30<tab>5565<tab>5671<tab>6928<tab>6821
    AE000111<tab>AE000394<tab>1e-27<tab>5587<tab>5671<tab>135<tab>219
    AE000111<tab>AE000425<tab>6e-26<tab>5587<tab>5671<tab>8552<tab>8468
    AE000111<tab>AE000171<tab>3e-24<tab>5587<tab>5671<tab>2214<tab>2130

An example Legacy BLAST 9-formatted file comparing two protein sequences, taken
from [3]_ (tab characters represented by ``<tab>``)::

    # BLASTN 2.2.3 [May-13-2002]
    # Database: other_vertebrate
    # Query: AF178033
    # Fields:
    Query id,Subject id,% identity,alignment length,mismatches,gap openings,q.\
 start,q. end,s. start,s. end,e-value,bit score
    AF178033<tab>EMORG:AF178033<tab>100.00<tab>811<tab>0<tab>0<tab>1<tab>811<t\
ab>1<tab>811<tab>0.0<tab>1566.6
    AF178033<tab>EMORG:AF031394<tab>99.63<tab>811<tab>3<tab>0<tab>1<tab>811<ta\
b>99<tab>909<tab>0.0<tab>1542.8
    AF178033<tab>EMORG:AF031393<tab>95.07<tab>811<tab>40<tab>0<tab>1<tab>811<t\
ab>99<tab>909<tab>0.0<tab>1249.4
    AF178033<tab>EMORG:AF178031<tab>94.82<tab>811<tab>42<tab>0<tab>1<tab>811<t\
ab>1<tab>811<tab>0.0<tab>1233.5
    AF178033<tab>EMORG:AF178032<tab>94.57<tab>811<tab>44<tab>0<tab>1<tab>811<t\
ab>1<tab>811<tab>0.0<tab>1217.7

Format Support
==============
**Has Sniffer: Yes**

**State: Experimental as of 0.4.0-dev.**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`pandas.DataFrame`                                        |
+------+------+---------------------------------------------------------------+

Format Specification
====================
Format Types
------------
BLAST+7
^^^^^^^
BLAST+7 format is a tabular text-based format produced by BLAST+ output
format 7 (``-outfmt 7``). Data is tab-separated and has no column headers.
Each comparison between a subject and either a query or a database is
displayed, even if there are no similarities. For example:

.. code-block:: none

    # BLASTP 2.2.31+
    # Query: query1
    # Subject: subject1
    # 0 hits found

Each pass through different combinations of query and subject contains the
first 3 lines shown above, but if similarities are found, the next line will
differ and other information will be added.
For example::

    # BLASTP 2.2.31+
    # Query: query1
    # Subject: subject2
    # Fields: q. start, q. end, s. start, s. end, identical, mismatches, sbjct\
frame, query acc.ver, subject acc.ver
    # 2 hits found
    1	8	3	10	8	0	1	query1	subject2
    2	5	2	15	8	0	2	query1	subject2

When data is found, a "Fields" line is added, expressing the types of returned
data, ordered respectively, the "hits found" line expresses how many rows of
data are returned, and at the bottom, the data is shown, separated by the tab
character.

Legacy BLAST 9
^^^^^^^^^^^^^^
Legacy BLAST 9 format is a tabular text-based format produced by Legacy
BLAST output format 9 (``-m 9``). Data is tab-separated and has no column
headers. Only one data block is shown, as ``blastall`` (the command that
returns Legacy BLAST data) only compares to a single database.

An example of a Legacy BLAST 9 file:

.. code-block:: none

    # BLASTN 2.2.3 [May-13-2002]
    # Database: other_vertebrate
    # Query: AF178033
    # Fields:
    Query id,Subject id,% identity,alignment length,mismatches,gap openings,q.\
 start,q. end,s. start,s. end,e-value,bit score
    AF178033    EMORG:AF178033  100.00  811 0   0   1   811 1   811 0.0 1566.6
    AF178033    EMORG:AF031394  99.63   811 3   0   1   811 99  909 0.0 1542.8

This format differs from BLAST+7 in a few ways: the fields are displayed on a
new line after "# Fields:", fields are not customizable, so they will always
be the default fields, Query id and Subject id are capitalized, and when
listing fields, there are no spaces after commas.

.. note:: Both data types ignore lines starting with "#", except for the
   "# Fields:" line in BLAST+7. They both only find the fields and data, as
   those are the only items present in the returned DataFrame.

.. note:: The reader supports both BLAST+7 and Legacy BLAST 9, and even a file
   with both present.

BLAST Column Types
------------------
The following column types are output by BLAST and supported by scikit-bio.
This table contains the outputs and their corresponding widely used name,
which is displayed on the returned DataFrame. For more information as to
meaning and data types of these fields, please see
:mod:`skbio.io.format.blast6`.

+-------------------+------------+
|Fields Name        |DataFrame C\|
|                   |olumn Name  |
+===================+============+
|query id           |qseqid      |
+-------------------+------------+
|query gi           |qgi         |
+-------------------+------------+
|query acc.         |qacc        |
+-------------------+------------+
|query acc.ver      |qaccver     |
+-------------------+------------+
|query length       |qlen        |
+-------------------+------------+
|subject id         |sseqid      |
+-------------------+------------+
|subject ids        |sallseqid   |
+-------------------+------------+
|subject gi         |sgi         |
+-------------------+------------+
|subject gis        |sallgi      |
+-------------------+------------+
|subject acc.       |sacc        |
+-------------------+------------+
|subject acc.ver    |saccver     |
+-------------------+------------+
|subject accs       |sallacc     |
+-------------------+------------+
|subject length     |slen        |
+-------------------+------------+
|q\. start          |qstart      |
+-------------------+------------+
|q\. end            |qend        |
+-------------------+------------+
|s\. start          |sstart      |
+-------------------+------------+
|s\. end            |send        |
+-------------------+------------+
|query seq          |qseq        |
+-------------------+------------+
|subject seq        |sseq        |
+-------------------+------------+
|evalue             |evalue      |
+-------------------+------------+
|bit score          |bitscore    |
+-------------------+------------+
|score              |score       |
+-------------------+------------+
|alignment length   |length      |
+-------------------+------------+
|% identity         |pident      |
+-------------------+------------+
|identical          |nident      |
+-------------------+------------+
|mismatches         |mismatch    |
+-------------------+------------+
|positives          |positive    |
+-------------------+------------+
|gap opens          |gapopen     |
+-------------------+------------+
|gaps               |gaps        |
+-------------------+------------+
|% positives        |ppos        |
+-------------------+------------+
|query/sbjct frames |frames      |
+-------------------+------------+
|query frame        |qframe      |
+-------------------+------------+
|sbjct frame        |sframe      |
+-------------------+------------+
|BTOP               |btop        |
+-------------------+------------+
|subject tax ids    |staxids     |
+-------------------+------------+
|subject sci names  |sscinames   |
+-------------------+------------+
|subject com names  |scomnames   |
+-------------------+------------+
|subject blast names|sblastnames |
+-------------------+------------+
|subject super king\|sskingdoms  |
|doms               |            |
+-------------------+------------+
|subject title      |stitle      |
+-------------------+------------+
|subject strand     |sstrand     |
+-------------------+------------+
|subject titles     |salltitles  |
+-------------------+------------+
|% query coverage p\|qcovs       |
|er subject         |            |
+-------------------+------------+
|% query coverage p\|qcovhsp     |
|er hsp             |            |
+-------------------+------------+

Examples
========
Suppose we have a BLAST+7 file:

>>> from io import StringIO
>>> import skbio.io
>>> import pandas as pd
>>> fs = '\\n'.join([
...     '# BLASTN 2.2.18+',
...     '# Query: gi|1786181|gb|AE000111.1|AE000111',
...     '# Database: ecoli',
...     '# Fields: query acc., subject acc., evalue, q. start, q. end, s. st\
art, s. end',
...     '# 5 hits found',
...     'AE000111\\tAE000111\\t0.0\\t1\\t10596\\t1\\t10596',
...     'AE000111\\tAE000174\\t8e-30\\t5565\\t5671\\t6928\\t6821',
...     'AE000111\\tAE000171\\t3e-24\\t5587\\t5671\\t2214\\t2130',
...     'AE000111\\tAE000425\\t6e-26\\t5587\\t5671\\t8552\\t8468'
... ])
>>> fh = StringIO(fs)

Read the file into a ``pd.DataFrame``:

>>> df = skbio.io.read(fh, format="blast+7", into=pd.DataFrame)
>>> df # doctest: +NORMALIZE_WHITESPACE
       qacc      sacc        evalue  qstart   qend  sstart   send
0  AE000111  AE000111  0.000000e+00       1  10596       1  10596
1  AE000111  AE000174  8.000000e-30    5565   5671    6928   6821
2  AE000111  AE000171  3.000000e-24    5587   5671    2214   2130
3  AE000111  AE000425  6.000000e-26    5587   5671    8552   8468

Suppose we have a Legacy BLAST 9 file:

>>> from io import StringIO
>>> import skbio.io
>>> import pandas as pd
>>> fs = '\\n'.join([
...     '# BLASTN 2.2.3 [May-13-2002]',
...     '# Database: other_vertebrate',
...     '# Query: AF178033',
...     '# Fields: ',
...     'Query id,Subject id,% identity,alignment length,mismatches,gap openin\
gs,q. start,q. end,s. start,s. end,e-value,bit score',
...     'AF178033\\tEMORG:AF178033\\t100.00\\t811\\t0\\t0\\t1\\t811\\t1\\t81\
1\\t0.0\\t1566.6',
...     'AF178033\\tEMORG:AF178032\\t94.57\\t811\\t44\\t0\\t1\\t811\\t1\\t81\
1\\t0.0\\t1217.7',
...     'AF178033\\tEMORG:AF178031\\t94.82\\t811\\t42\\t0\\t1\\t811\\t1\\t81\
1\\t0.0\\t1233.5'
... ])
>>> fh = StringIO(fs)

Read the file into a ``pd.DataFrame``:

>>> df = skbio.io.read(fh, format="blast+7", into=pd.DataFrame)
>>> df # doctest: +NORMALIZE_WHITESPACE
     qseqid          sseqid  pident  length  mismatch  gapopen  qstart  qend \\
0  AF178033  EMORG:AF178033  100.00     811         0        0       1   811
1  AF178033  EMORG:AF178032   94.57     811        44        0       1   811
2  AF178033  EMORG:AF178031   94.82     811        42        0       1   811
<BLANKLINE>
   sstart  send  evalue  bitscore
0       1   811       0    1566.6
1       1   811       0    1217.7
2       1   811       0    1233.5

References
==========
.. [1] Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990)
   "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
.. [2] http://www.ncbi.nlm.nih.gov/books/NBK279682/
.. [3] http://www.compbio.ox.ac.uk/analysis_tools/BLAST/BLAST_blastall/blastal\
l_examples.shtml
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
import numpy as np

from skbio.io import create_format, BLAST7FormatError
from skbio.io.format.blast6 import _possible_columns

blast7 = create_format('blast+7')

column_converter = {'query id': 'qseqid', 'query gi': 'qgi',
                    'query acc.': 'qacc', 'query acc.ver': 'qaccver',
                    'query length': 'qlen', 'subject id': 'sseqid',
                    'subject ids': 'sallseqid', 'subject gi': 'sgi',
                    'subject gis': 'sallgi', 'subject acc.': 'sacc',
                    'subject acc.ver': 'saccver', 'subject accs.': 'sallacc',
                    'subject length': 'slen', 'q. start': 'qstart',
                    'q. end': 'qend', 's. start': 'sstart', 's. end': 'send',
                    'query seq': 'qseq', 'subject seq': 'sseq',
                    'evalue': 'evalue', 'bit score': 'bitscore',
                    'score': 'score', 'alignment length': 'length',
                    '% identity': 'pident', 'identical': 'nident',
                    'mismatches': 'mismatch', 'positives': 'positive',
                    'gap opens': 'gapopen', 'gaps': 'gaps',
                    '% positives': 'ppos', 'query/sbjct frames': 'frames',
                    'query frame': 'qframe', 'sbjct frame': 'sframe',
                    'BTOP': 'btop', 'subject tax ids': 'staxids',
                    'subject sci names': 'sscinames',
                    'subject com names': 'scomnames',
                    'subject blast names': 'sblastnames',
                    'subject super kingdoms': 'sskingdoms',
                    'subject title': 'stitle', 'subject titles': 'salltitles',
                    'subject strand': 'sstrand',
                    '% query coverage per subject': 'qcovs',
                    '% query coverage per hsp': 'qcovhsp',
                    'Query id': 'qseqid', 'Subject id': 'sseqid',
                    'gap openings': 'gapopen', 'e-value': 'evalue'}


@blast7.sniffer()
def _blast7_sniffer(fh):
    # Smells a BLAST+7 file if the following conditions are present
    #   -First line contains "BLAST"
    #   -Second line contains "Query" or "Database"
    #   -Third line starts with "Subject" or "Query"
    lines = list(zip(range(3), fh))
    if len(lines) < 3:
        return False, {}
    try:
        if (lines[0][1][:7] == "# BLAST" and
            (lines[1][1][:8] == "# Query:" or
             lines[1][1][:11] == "# Database:") and
            (lines[2][1][:10] == "# Subject:" or
             lines[2][1][:8] == "# Query:" or
             lines[2][1][:11] == "# Database:")):
            return True, {}
        else:
            return False, {}
    except IndexError:
        return False, {}


@blast7.reader(pd.DataFrame, monkey_patch=False)
def _blast7_to_data_frame(fh):
    dtypes = None
    x = 0
    y = 0
    skiprows = []
    for line in fh:
        if x == 1 and not dtypes:
            dtypes = _remove_newline_give_dtypes(line, x=x)
            x = 0
            skiprows.append(y)
        elif x == 1:
            next_dtypes = _remove_newline_give_dtypes(line, x=x)
            if dtypes != next_dtypes:
                raise BLAST7FormatError("Fields %r do not equal fields %r"
                                        % (dtypes, next_dtypes))
            x = 0
            skiprows.append(y)
        if line == "# Fields: \n":
            # Identifies Legacy BLAST 9 data
            x = 1
        elif line[:10] == "# Fields: ":
            # Identifies BLAST+7 data
            x = 2
        if x == 2 and not dtypes:
            dtypes = _remove_newline_give_dtypes(line)
            x = 0
        elif x == 2:
            # Affirms data types do not differ throught file
            next_dtypes = _remove_newline_give_dtypes(line)
            if dtypes != next_dtypes:
                raise BLAST7FormatError("Fields %r do not equal fields %r"
                                        % (dtypes, next_dtypes))
            x = 0
        y += 1
    if not dtypes:
        # Affirms file contains BLAST data
        raise BLAST7FormatError("File contains no BLAST data.")
    fh.seek(0)
    try:
        df = pd.read_csv(fh, names=dtypes, sep='\t', comment='#',
                         dtype=_possible_columns, na_values='N/A',
                         keep_default_na=False, skiprows=skiprows)
    except ValueError:
        # Catches first possibility of incorrect number of columns
        raise BLAST7FormatError("Number of data columns does not equal number"
                                " of fields (%r)." % (len(dtypes)))
    if df.index.dtype != np.int64:
        # Catches second possibility of incorrect number of columns
        raise BLAST7FormatError("Number of data columns does not equal number"
                                " of fields (%r)." % (len(dtypes)))
    return df


def _remove_newline_give_dtypes(line, x=2):
    # Removes "\n" from fields line and returns fields as a list (dtypes)
    line = line.rstrip('\n')
    if x == 1:
        fields = line.split(',')
    elif x == 2:
        line = line[10:]
        fields = line.split(', ')
    try:
        dtypes = [column_converter[field] for field in fields]
    except KeyError as e:
        # Affirms all column names are supported
        raise BLAST7FormatError("Unrecognized field: %r."
                                " Supported fields: %r"
                                % (str(e).strip("'"),
                                   set(_possible_columns.keys())))
    return dtypes
