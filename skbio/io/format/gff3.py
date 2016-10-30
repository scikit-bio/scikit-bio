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
|Yes   |Yes   |:mod:`skbio.metadata.IntervalMetadata`                         |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |generator of :mod:`skbio.metadata.IntervalMetadata`            |
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

The columns and attributes are read in as the vocabulary defined in
genbank parsers (:mod:`skbio.io.format.genbank`).

Examples
--------

>>> gff_str = """
... ##gff-version\t3.2.1
... ##sequence-region\tctg123\t1\t1497228
... ctg123\t.\tgene\t1000\t9000\t.\t+\t0\tID=gene00001;Name=EDEN
... ctg123\t.\tTF_binding_site\t1000\t1012\t.\t+\t.\tParent=gene00001
... ctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001;Parent=gene00001
... """
>>> import io
>>> from skbio.metadata import IntervalMetadata
>>> from skbio.io import read
>>> im = IntervalMetadata(20000)
>>> gff = io.StringIO(gff_str)
>>> im_return = read(gff, format='gff3', into=IntervalMetadata,
...                  interval_metadata=im)
>>> im_return   # doctest: +SKIP
3 interval features
-------------------
Interval(interval_metadata=<4601272528>, bounds=[(999, 9000)], fuzzy=\
[(False, False)], metadata={'source': '.', 'type': 'gene', 'strand': '+', \
'score': '.', 'phase': 0, 'ID': 'gene00001', 'Name': 'EDEN'})
Interval(interval_metadata=<4601272528>, bounds=[(999, 1012)], fuzzy=\
[(False, False)], metadata={'source': '.', 'type': 'TF_binding_site', \
'strand': '+', 'score': '.', 'Parent': 'gene00001'})
Interval(interval_metadata=<4601272528>, bounds=[(1049, 9000)], fuzzy=\
[(False, False)], metadata={'source': '.', 'type': 'mRNA', 'strand': '+', \
'score': '.', 'ID': 'mRNA00001', 'Parent': 'gene00001'})
>>> im == im_return
True

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
from io import StringIO

from skbio.sequence import DNA, Sequence
from skbio.io import create_format, GFF3FormatError
from skbio.metadata import IntervalMetadata
from skbio.io.format._base import (
    _line_generator, _too_many_blanks)
from skbio.io.format._base import _get_nth_sequence as _get_nth_record
from skbio.io.format._sequence_feature_vocabulary import (
    _vocabulary_change, _vocabulary_skip)
from skbio.io import write, read


gff3 = create_format('gff3')


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
def _gff3_to_generator(fh, interval_metadata_dict):
    '''Parse the GFF3 into the existing IntervalMetadata

    Note that if the seq ID does not exist in the input dict of
    `interval_metadata_dict`, the record will be skipped and not
    parsed.

    Parameters
    ----------
    fh : file
        file handler
    interval_metadata_dict : dict
        key is seq ID and value is the IntervalMetadata for the seq.

    Yields
    ------
    IntervalMetadata
    '''
    for seq_id, lines in _yield_record(fh):
        if seq_id in interval_metadata_dict:
            imd = interval_metadata_dict[seq_id]
            yield _parse_record(lines, imd)


@gff3.writer(None)
def _generator_to_gff3(obj, fh, seq_ids, skip=True):
    '''Write list of IntervalMetadata into file.

    Parameters
    ----------
    obj : Iterable of IntervalMetadata
    fh : file handler
    seq_ids : Iterable of seq id (str)
    '''
    for obj_i, seq_id in zip(obj, seq_ids):
        _serialize_interval_metadata(obj_i, seq_id, fh, skip)


@gff3.reader(Sequence)
def _gff3_to_sequence(fh, rec_num=1):
    ''''''
    seq_id, lines = _get_nth_record(_yield_record(fh), rec_num)
    # you can't read directly from fh because in registry.py line 543
    # file.tell() will fail "telling position disabled by next() call".
    stream = StringIO(fh.read())
    seq = read(stream, format='fasta', into=Sequence, seq_num=rec_num)
    _parse_record(lines, interval_metadata=seq.interval_metadata)
    return seq


@gff3.writer(Sequence)
def _sequence_to_gff3(obj, fh, skip=True):
    _serialize_seq(obj, fh, skip)


@gff3.reader(DNA)
def _gff3_to_dna(fh, rec_num=1):
    ''''''
    seq_id, lines = _get_nth_record(_yield_record(fh), rec_num)
    stream = StringIO(fh.read())
    seq = read(stream, format='fasta', into=DNA, seq_num=rec_num)
    _parse_record(lines, interval_metadata=seq.interval_metadata)
    return seq


@gff3.writer(DNA)
def _dna_to_gff3(obj, fh, skip=True):
    _serialize_seq(obj, fh, skip)


@gff3.reader(IntervalMetadata)
def _gff3_to_interval_metadata(fh, interval_metadata, rec_num=1):
    '''Read a GFF3 record into the specified interval metadata.

    Parameters
    ----------
    fh : file handler
    interval_metadata : IntervalMetadata
    rec_num : int
        which record to read in.
    '''
    seq_id, lines = _get_nth_record(_yield_record(fh), rec_num)
    return _parse_record(lines, interval_metadata=interval_metadata)


@gff3.writer(IntervalMetadata)
def _interval_metadata_to_gff3(obj, fh, seq_id, skip=True):
    '''
    Parameters
    ----------
    obj : IntervalMetadata
    seq_id : str
        ID for column 1 in the GFF3 file.
    '''
    _serialize_interval_metadata(obj, seq_id, fh, skip=True)


def _yield_record(fh):
    '''Yield lines that belong to the same sequence.'''
    lines = []
    current = False
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        if line.startswith('##FASTA'):
            # stop once reaching to sequence section
            break
        if not line.startswith('#'):
            try:
                seq_id, _ = line.split('\t', 1)
            except ValueError:
                raise GFF3FormatError(
                    'Wrong GFF3 format at line: %s' % line)
            if current == seq_id:
                lines.append(line)
            else:
                if current is not False:
                    yield current, lines
                lines = [line]
                current = seq_id
    yield current, lines


def _parse_record(lines, interval_metadata):
    '''Parse the lines into a IntervalMetadata object.'''
    for line in lines:
        columns = line.split('\t')
        # there should be 9 columns
        if len(columns) != 9:
            raise GFF3FormatError(
                'do not have 9 columns in this line: "%s"' % line)
        # the 1st column is seq ID for every feature. don't store
        # this repetitive information
        metadata = {'source': columns[1],
                    'type': columns[2],
                    'score': columns[5],
                    'strand': columns[6]}
        phase = columns[7]
        try:
            metadata['phase'] = int(phase)
        except ValueError:
            if phase != '.':
                raise GFF3FormatError(
                    'unknown value for phase column: {!r}'.format(phase))
        metadata.update(_parse_attr(columns[8]))

        start, end = columns[3:5]

        bounds = [(int(start)-1, int(end))]

        interval_metadata.add(bounds, metadata=metadata)

    return interval_metadata


def _parse_attr(s):
    '''parse attribute column'''
    voca_change = _vocabulary_change('gff3')
    md = {}
    # in case the line ending with ';', strip it.
    s = s.rstrip(';')
    for attr in s.split(';'):
        k, v = attr.split('=')
        if k in voca_change:
            k = voca_change[k]
        md[k] = v
    return md


def _serialize_interval_metadata(
        interval_metadata, seq_id, fh, skip=True):
    '''Serialize an IntervalMetadata to GFF3.

    Parameters
    ----------
    interval_metadata : IntervalMetadata
    skip : bool
        whether to skip outputting each sub region as a line in GFF3.
    '''
    # write file header
    print('##gff-version 3', file=fh)

    column_keys = ['source', 'type', 'score', 'strand', 'phase']
    voca_change = _vocabulary_change('gff3', False)
    voca_skip = _vocabulary_skip('gff3')
    voca_skip.extend(column_keys)

    # these characters have reserved meanings in column 9 and must be
    # escaped when used in other contexts
    escape = str.maketrans({';': '%3B',
                            '=': '%3D',
                            '&': '%26',
                            ',': '%2C'})

    for interval in interval_metadata._intervals:
        md = interval.metadata
        bd = interval.bounds
        start = str(bd[0][0] + 1)
        end = str(bd[-1][-1])

        source, feat_type, score, strand, phase = [
            str(md.get(i, '.')) for i in column_keys]
        columns = [seq_id, source, feat_type, start, end, score, strand, phase]

        # serialize the attributes in column 9
        attr = []
        # use sort to make the output order deterministic
        for k in sorted(md):
            if k in voca_skip:
                # skip the metadata don't go to attribute column
                continue
            if k in voca_change:
                k = voca_change[k]
            v = md[k]
            attr.append('%s=%s' % (k.translate(escape), v.translate(escape)))
        columns.append(';'.join(attr))

        print('\t'.join(columns), file=fh)

        # if there are multiple regions for this feature,
        # output each region as a standalone line in GFF3.
        if len(bd) > 1 and skip is False:
            for start, end in bd:
                # if this is a gene, then each sub region should be an exon
                if columns[2] == 'gene':
                    columns[2] = 'exon'
                columns[3] = str(start + 1)
                columns[4] = str(end)
                try:
                    parent = md['ID']
                except KeyError:
                    raise GFF3FormatError(
                        'You need provide ID info for '
                        'the parent interval feature: %r' % interval)
                columns[8] = 'Parent=%s' % parent
                print('\t'.join(columns), file=fh)


def _serialize_seq(seq, fh, skip=True):
    '''Serialize a sequence to GFF3.'''
    _serialize_interval_metadata(
        seq.interval_metadata, seq.metadata['id'], fh, skip)
    print('##FASTA', file=fh)
    write(seq, into=fh, format='fasta')
