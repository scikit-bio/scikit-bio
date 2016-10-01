r'''GFF3 format (:mod:`skbio.io.format.gff3)
============================================

GFF3 is a standard file format for storing genomic features in a text
file. GFF stands for Generic Feature Format. GFF files are plain
text, 9 column, tab-delimited files [#]_.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.metadata.IntervalMetadata` objects                 |
+------+------+---------------------------------------------------------------+
|Yes   |No    |generator of :mod:`skbio.metadata.IntervalMetadata` objects    |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
The first line of the file is a comment that identifies the format and
version.  This is followed by a series of data lines, each one of
which corresponds to an annotation.  The 9 columns of the annotation
section are as follows:

+--------+------------------------------------------+
| Column | Description                              |
+========+==========================================+
| SEQID  | ID of the landmark used                  |
+--------+------------------------------------------+
| SOURCE | algorithm used to generate this feature  |
+--------+------------------------------------------+
| TYPE   | type of the feature                      |
+--------+------------------------------------------+
| START  | start of the feature                     |
+--------+------------------------------------------+
| END    | end of the feature                       |
+--------+------------------------------------------+
| SCORE  | floating point score                     |
+--------+------------------------------------------+
| STRAND | The strand of the feature (+/-/./?)      |
+--------+------------------------------------------+
| PHASE  | only for TYPE="CDS"                      |
+--------+------------------------------------------+
| ATTR   | feature attributes                       |
+--------+------------------------------------------+


Column 9 (ATTR) is list of feature attributes in the format
tag=value. Multiple tag=value pairs are separated by
semicolons. Multiple values of the same tag are indicated by
separating the values with the comma ",". The following tags have
predefined meanings:

* ID. Indicates the unique identifier of the feature. IDs must be
  unique within the scope of the GFF file.

* Name. Display name for the feature. This is the name to be displayed
  to the user. Unlike IDs, there is no requirement that the Name be
  unique within the file.

* Alias. A secondary name for the feature. It is suggested that this
  tag be used whenever a secondary identifier for the feature is
  needed, such as locus names and accession numbers. Unlike ID, there
  is no requirement that Alias be unique within the file.

* Parent. Indicates the parent of the feature. A parent ID can be used
  to group exons into transcripts, transcripts into genes and so
  forth. A feature may have multiple parents. Parent can *only* be
  used to indicate a partof relationship.

* Target. Indicates the target of a nucleotide-to-nucleotide or
  protein-to-nucleotide alignment. The format of the value is
  "target_id start end [strand]", where strand is optional and may be
  "+" or "-". If the target_id contains spaces, they must be escaped
  as hex escape %20.

* Gap. The alignment of the feature to the target if the two are not
  collinear (e.g. contain gaps). The alignment format is taken from
  the CIGAR format described in the Exonerate documentation.

* Derives_from. Used to disambiguate the relationship between one
  feature and another when the relationship is a temporal one rather
  than a purely structural "part of" one. This is needed for
  polycistronic genes.

* Note. A free text note.

* Dbxref. A database cross reference. See the GFF3 specification for
  more information.

* Ontology_term. A cross reference to an ontology term.

* Is_circular. A flag to indicate whether a feature is circular.

Examples
--------

>>> gff_str = """
... ##gff-version\t3.2.1
... ##sequence-region\tctg123\t1\t1497228
... ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN
... ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tParent=gene00001
... ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001
... """
>>> import io
>>> from skbio.metadata import IntervalMetadata
>>> from skbio.io import read
>>> gff = io.StringIO(gff_str)
>>> imd = read(gff, format='gff3', into=IntervalMetadata, upper_bound=2000000)
>>> imd   # doctest: +SKIP
3 interval features
-------------------
Interval(interval_metadata=<4601272528>, bounds=[(999, 9000)], fuzzy=\
[(False, False)], metadata={'SOURCE': '.', 'TYPE': 'gene', 'STRAND': '+', \
'SCORE': '.', 'PHASE': '.', 'ATTR': 'ID=gene00001;Name=EDEN'})
Interval(interval_metadata=<4601272528>, bounds=[(999, 1012)], fuzzy=\
[(False, False)], metadata={'SOURCE': '.', 'TYPE': 'TF_binding_site', \
'STRAND': '+', 'SCORE': '.', 'PHASE': '.', 'ATTR': 'Parent=gene00001'})
Interval(interval_metadata=<4601272528>, bounds=[(1049, 9000)], fuzzy=\
[(False, False)], metadata={'SOURCE': '.', 'TYPE': 'mRNA', 'STRAND': '+', \
'SCORE': '.', 'PHASE': '.', 'ATTR': 'ID=mRNA00001;Parent=gene00001'})

Reference
---------

.. [#] https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md  # noqa

'''

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re

from skbio.io import create_format, FileFormatError
from skbio.metadata import IntervalMetadata
from skbio.io.format._base import (
    _line_generator, _too_many_blanks)
from skbio.io.format._base import _get_nth_sequence as _get_nth_record


class GFF3FormatError(FileFormatError):
    pass


gff3 = create_format('gff3')


_GFF3_HEADERS = [
    'SEQID',
    'SOURCE',
    'TYPE',
    'START',
    'END',
    'SCORE',
    'STRAND',
    'PHASE',
    'ATTR']


@gff3.sniffer()
def _gff3_sniffer(fh):
    # check the 1st real line is a valid ID line
    if _too_many_blanks(fh, 5):
        return False, {}

    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    if re.match(r'##gff-version\s+3', line):
        return True, {}
    else:
        return False, {}


@gff3.reader(None)
def _gff3_to_generator(fh, upper_bounds):
    return _parse_records(fh, upper_bounds)


@gff3.reader(IntervalMetadata)
def _gff3_to_interval_metadata(fh, upper_bound, rec_num=1):
    return _get_nth_record(_parse_records(fh, [upper_bound]), rec_num)


def _parse_records(fh, upper_bounds):
    '''Yield record

    A record is all the lines associated with the same ID (seq). This
    parser does minimal parsing. It does not convert variable types or
    parse the attr of column 9.

    '''
    seq_id = False
    i = 0
    imd = IntervalMetadata(upper_bounds[i])
    md_headers = _GFF3_HEADERS[1:3] + _GFF3_HEADERS[5:]
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        if not line.startswith('#'):
            columns = line.split('\t')
            if len(columns) != len(_GFF3_HEADERS):
                raise GFF3FormatError(
                    'do not have 9 columns in this line: "%s"' % line)
            # the 1st column is seq ID for every feature. don't store
            # this repetitive information
            metadata = dict(zip(md_headers,
                                columns[1:3] + columns[5:]))

            start, end = columns[3:5]

            bounds = [(int(start)-1, int(end))]

            if columns[0] == seq_id or seq_id is False:
                seq_id = columns[0]
                imd.add(bounds, metadata=metadata)
            else:
                seq_id = False
                yield imd
                i += 1
                imd = IntervalMetadata(upper_bounds[i])
                imd.add(bounds, metadata=metadata)

    yield imd


@gff3.writer(IntervalMetadata)
def _interval_metadata_to_gff3(obj, fh, seq_id):
    '''
    Parameters
    ----------
    obj : IntervalMetadata
    seq_id : str
        ID for column 1 in the GFF3 file.
    '''
    # write file header
    print('##gff-version 3', file=fh)

    for interval in obj._intervals:
        columns = [seq_id]
        md = interval.metadata
        for bound in interval.bounds:
            for head in _GFF3_HEADERS[1:]:
                if head == 'START':
                    columns.append(str(bound[0]+1))
                elif head == 'END':
                    columns.append(str(bound[1]))
                elif head in md:
                    columns.append(md[head])
                else:
                    columns.append('.')
            print('\t'.join(columns), file=fh)
