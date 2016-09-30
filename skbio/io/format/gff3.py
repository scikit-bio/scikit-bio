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
... ##gff-version	3.2.1
... ##sequence-region	ctg123	1	1497228
... ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN
... ctg123	.	TF_binding_site	1000	1012	.	+	.	Parent=gene00001
... ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001
... """
>>> import io
>>> from skbio.metadata import IntervalMetadata
>>> from skbio.io import read
>>> gff = io.StringIO(gff_str)
>>> imd = read(gff, format='gff3', into=IntervalMetadata, upper_bound=2000000)
>>> imd

Reference
---------

.. [#] https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

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

    if re.match(r'##gff-version +3', line):
        return True, {}
    else:
        return False, {}


def _is_float(input):
    try:
        float(input)
    except ValueError:
        return False
    return True


def _parse_required(s):
    if s.isdigit():
        return int(s)
    elif _is_float(s):
        return float(s)
    else:
        return s


@gff3.reader(None)
def _gff3_to_generator(fh, upper_bounds):
    return _parse_records(fh, upper_bounds)


@gff3.reader(IntervalMetadata)
def _gff3_to_interval_metadata(fh, upper_bound, rec_num=1):
    return  _get_nth_record(_parse_records(fh, [upper_bound]), rec_num)


def _parse_records(fh, upper_bounds):
    '''Yield record

    A record is all the lines associated with the same ID (seq).
    '''
    seq_id = False
    i = 0
    imd = IntervalMetadata(upper_bounds[i])
    md_headers = _GFF3_HEADERS[1:3] + _GFF3_HEADERS[5:]
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        if not line.startswith('#'):
            columns = line.split('\t')
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


def _attr_to_list(attr_list):
    _tags = []
    if any(isinstance(el, tuple) for el in attr_list):
        for _attr in attr_list:
            _tags.append('='.join([_attr[0], _attr[1]]))
        return ';'.join(_tags)
    else:
        return '='.join([attr_list[0], attr_list[1]])


@gff3.writer(IntervalMetadata)
def _interval_metadata_to_gff3(obj, fh):
    # write file header
    fh.write('%s\n' % '##gff-version 3')

    for features, interval in obj.features.items():
        tab = [None] * 9

        if len(features) is not 7:
            raise GFF3FormatError(
                "``IntervalMetadata`` can only be written in GFF3 format if all"
                " annotation columns are found.")
        if len(interval) is not 2:
            raise GFF3FormatError(
                "``IntervalMetadata`` can only be written in GFF3 format if "
                " `START` and `END` fields are provided.")
        if not all(_annot in set(features) for _annot in _ANNOTATION_HEADERS):
            raise GFF3FormatError(
                "GFF3 format requires header names to match pre-defined set: %s"
                % ', '.join(_ANNOTATION_HEADERS))

        for i, annot in enumerate(_GFF3_HEADERS):
            if annot is 'START':
                tab[i] = interval[0]
            elif annot is 'END':
                tab[i] = interval[1]
            elif annot is 'ATTR':
                tab[i] = _attr_to_list(features[annot])
            else:
                tab[i] = features[annot]

        fh.write('\t'.join(map(str, tab)) + '\n')
