"""
BLAST+7 format (:mod:`skbio.io.format.blast7`)
==============================================

.. currentmodule:: skbio.io.format.blast7

The BLAST+7 format (``blast+7``) stores the results of a BLAST [1]_ database
search. This format is produced by both BLAST+ output
format 7 and legacy BLAST output format 9. The results
are stored in a simple tabular format with headers. Values are separated by the
tab character.

An example BLAST+7-formatted file comparing two nucleotide sequences, taken
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

Format Support
==============
**Has Sniffer: Yes**

**State: Experimental as of 0.4.1.**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |:mod:`pandas.DataFrame`                                        |
+------+------+---------------------------------------------------------------+

Format Specification
====================
There are two BLAST+7 file formats supported by scikit-bio: BLAST+ output
format 7 (``-outfmt 7``) and legacy BLAST output format 9 (``-m 9``). Both file
formats are structurally similar, with minor differences.

Example BLAST+ output format 7 file::

    # BLASTP 2.2.31+
    # Query: query1
    # Subject: subject2
    # Fields: q. start, q. end, s. start, s. end, identical, mismatches, sbjct\
frame, query acc.ver, subject acc.ver
    # 2 hits found
    1	8	3	10	8	0	1	query1	subject2
    2	5	2	15	8	0	2	query1	subject2

.. note:: Database searches without hits may occur in BLAST+ output format 7
   files. scikit-bio ignores these "empty" records:

   .. code-block:: none

       # BLASTP 2.2.31+
       # Query: query1
       # Subject: subject1
       # 0 hits found

Example legacy BLAST output format 9 file:

.. code-block:: none

    # BLASTN 2.2.3 [May-13-2002]
    # Database: other_vertebrate
    # Query: AF178033
    # Fields:
    Query id,Subject id,% identity,alignment length,mismatches,gap openings,q.\
 start,q. end,s. start,s. end,e-value,bit score
    AF178033    EMORG:AF178033  100.00  811 0   0   1   811 1   811 0.0 1566.6
    AF178033    EMORG:AF031394  99.63   811 3   0   1   811 99  909 0.0 1542.8

.. note:: scikit-bio requires fields to be consistent within a file.

BLAST Column Types
------------------
The following column types are output by BLAST and supported by scikit-bio.
For more information on these column types, see :mod:`skbio.io.format.blast6`.

+-------------------+----------------------+
|Field Name         |DataFrame Column Name |
+===================+======================+
|query id           |qseqid                |
+-------------------+----------------------+
|query gi           |qgi                   |
+-------------------+----------------------+
|query acc.         |qacc                  |
+-------------------+----------------------+
|query acc.ver      |qaccver               |
+-------------------+----------------------+
|query length       |qlen                  |
+-------------------+----------------------+
|subject id         |sseqid                |
+-------------------+----------------------+
|subject ids        |sallseqid             |
+-------------------+----------------------+
|subject gi         |sgi                   |
+-------------------+----------------------+
|subject gis        |sallgi                |
+-------------------+----------------------+
|subject acc.       |sacc                  |
+-------------------+----------------------+
|subject acc.ver    |saccver               |
+-------------------+----------------------+
|subject accs       |sallacc               |
+-------------------+----------------------+
|subject length     |slen                  |
+-------------------+----------------------+
|q\\. start          |qstart                |
+-------------------+----------------------+
|q\\. end            |qend                  |
+-------------------+----------------------+
|s\\. start          |sstart                |
+-------------------+----------------------+
|s\\. end            |send                  |
+-------------------+----------------------+
|query seq          |qseq                  |
+-------------------+----------------------+
|subject seq        |sseq                  |
+-------------------+----------------------+
|evalue             |evalue                |
+-------------------+----------------------+
|bit score          |bitscore              |
+-------------------+----------------------+
|score              |score                 |
+-------------------+----------------------+
|alignment length   |length                |
+-------------------+----------------------+
|% identity         |pident                |
+-------------------+----------------------+
|identical          |nident                |
+-------------------+----------------------+
|mismatches         |mismatch              |
+-------------------+----------------------+
|positives          |positive              |
+-------------------+----------------------+
|gap opens          |gapopen               |
+-------------------+----------------------+
|gaps               |gaps                  |
+-------------------+----------------------+
|% positives        |ppos                  |
+-------------------+----------------------+
|query/sbjct frames |frames                |
+-------------------+----------------------+
|query frame        |qframe                |
+-------------------+----------------------+
|sbjct frame        |sframe                |
+-------------------+----------------------+
|BTOP               |btop                  |
+-------------------+----------------------+
|subject tax ids    |staxids               |
+-------------------+----------------------+
|subject sci names  |sscinames             |
+-------------------+----------------------+
|subject com names  |scomnames             |
+-------------------+----------------------+
|subject blast names|sblastnames           |
+-------------------+----------------------+
|subject super      |sskingdoms            |
|kingdoms           |                      |
+-------------------+----------------------+
|subject title      |stitle                |
+-------------------+----------------------+
|subject strand     |sstrand               |
+-------------------+----------------------+
|subject titles     |salltitles            |
+-------------------+----------------------+
|% query coverage   |qcovs                 |
|per subject        |                      |
+-------------------+----------------------+
|% query coverage   |qcovhsp               |
|per hsp            |                      |
+-------------------+----------------------+

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

>>> df = skbio.io.read(fh, into=pd.DataFrame)
>>> df # doctest: +NORMALIZE_WHITESPACE
       qacc      sacc        evalue  qstart     qend  sstart     send
0  AE000111  AE000111  0.000000e+00     1.0  10596.0     1.0  10596.0
1  AE000111  AE000174  8.000000e-30  5565.0   5671.0  6928.0   6821.0
2  AE000111  AE000171  3.000000e-24  5587.0   5671.0  2214.0   2130.0
3  AE000111  AE000425  6.000000e-26  5587.0   5671.0  8552.0   8468.0

Suppose we have a legacy BLAST 9 file:

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

>>> df = skbio.io.read(fh, into=pd.DataFrame)
>>> df # doctest: +NORMALIZE_WHITESPACE
     qseqid          sseqid  pident  length  mismatch  gapopen  qstart  qend \\
0  AF178033  EMORG:AF178033  100.00   811.0       0.0      0.0     1.0  811.0
1  AF178033  EMORG:AF178032   94.57   811.0      44.0      0.0     1.0  811.0
2  AF178033  EMORG:AF178031   94.82   811.0      42.0      0.0     1.0  811.0
<BLANKLINE>
   sstart   send  evalue  bitscore
0     1.0  811.0     0.0    1566.6
1     1.0  811.0     0.0    1217.7
2     1.0  811.0     0.0    1233.5

References
==========
.. [1] Altschul, S.F., Gish, W., Miller, W., Myers, E.W. & Lipman, D.J. (1990)
   "Basic local alignment search tool." J. Mol. Biol. 215:403-410.
.. [2] http://www.ncbi.nlm.nih.gov/books/NBK279682/
"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

from skbio.io import create_format, BLAST7FormatError
from skbio.io.format._blast import _parse_blast_data

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
    #   -Third line starts with "Subject" or "Query" or "Database"
    lines = [line for _, line in zip(range(3), fh)]
    if len(lines) < 3:
        return False, {}

    if not lines[0].startswith("# BLAST"):
        return False, {}
    if not (lines[1].startswith("# Query:") or
            lines[1].startswith("# Database:")):
        return False, {}
    if not (lines[2].startswith("# Subject:") or
            lines[2].startswith("# Query:") or
            lines[2].startswith("# Database:")):
        return False, {}

    return True, {}


@blast7.reader(pd.DataFrame, monkey_patch=False)
def _blast7_to_data_frame(fh):
    line_num = 0
    columns = None
    skiprows = []
    for line in fh:
        if line == "# Fields: \n":
            # Identifies Legacy BLAST 9 data
            line = next(fh)
            line_num += 1
            if columns is None:
                columns = _parse_fields(line, legacy=True)
                skiprows.append(line_num)
            else:
                next_columns = _parse_fields(line, legacy=True)
                if columns != next_columns:
                    raise BLAST7FormatError("Fields %r do not equal fields %r"
                                            % (columns, next_columns))
                skiprows.append(line_num)
        elif line.startswith("# Fields: "):
            # Identifies BLAST+7 data
            if columns is None:
                columns = _parse_fields(line)
            else:
                # Affirms data types do not differ throught file
                next_columns = _parse_fields(line)
                if columns != next_columns:
                    raise BLAST7FormatError("Fields %r do not equal fields %r"
                                            % (columns, next_columns))
        line_num += 1
    if columns is None:
        # Affirms file contains BLAST data
        raise BLAST7FormatError("File contains no BLAST data.")
    fh.seek(0)

    return _parse_blast_data(fh, columns, BLAST7FormatError,
                             "Number of fields (%r) does not equal number"
                             " of data columns (%r).", comment='#',
                             skiprows=skiprows)


def _parse_fields(line, legacy=False):
    """Removes '\n' from fields line and returns fields as a list (columns)."""
    line = line.rstrip('\n')
    if legacy:
        fields = line.split(',')
    else:
        line = line.split('# Fields: ')[1]
        fields = line.split(', ')
    columns = []
    for field in fields:
        if field not in column_converter:
            raise BLAST7FormatError("Unrecognized field (%r)."
                                    " Supported fields: %r"
                                    % (field,
                                       set(column_converter.keys())))
        columns.append(column_converter[field])
    return columns
