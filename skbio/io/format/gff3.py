'''
GFF3 format (:mod:`skbio.io.format.gff3`)
=========================================

.. currentmodule:: skbio.io.format.gff3

GFF3 is a standard file format for storing genomic features in a text
file. GFF stands for Generic Feature Format. GFF files are plain
text, 9-column, tab-delimited files [1]_.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.metadata.IntervalMetadata`                         |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |generator of tuple (seq_id of str type,                        |
|      |      |:mod:`skbio.metadata.IntervalMetadata`)                        |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
**State: Experimental as of 0.5.1.**

The first line of the file is a comment that identifies the format and
version. This is followed by a series of data lines. Each data line
corresponds to an annotation and consists of 9 columns: SEQID, SOURCE,
TYPE, START, END, SCORE, STRAND, PHASE, and ATTR.

Column 9 (ATTR) is list of feature attributes in the format
"tag=value". Multiple "tag=value" pairs are delimited by
semicolons. Multiple values of the same tag are separated with the
comma ",". The following tags have predefined meanings: ID, Name,
Alias, Parent, Target, Gap, Derives_from, Note, Dbxref, Ontology_term,
and Is_circular.

The meaning and format of these columns and attributes are explained
detail in the format specification [1]_. And they are read in as the
vocabulary defined in GenBank parser (:mod:`skbio.io.format.genbank`).

Format Parameters
-----------------

Reader-specific Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^
Reader of ``IntervalMetadata`` generator requires a parameter
``lengths``. It is an iterable of int and specifies the upper bounds
of the ``IntervalMetadata`` objects the reader yields.

``IntervalMetadata`` GFF3 reader requires 2 parameters: ``seq_id`` of
str and ``length`` of int. It reads the annotation with the specified
sequence ID from the GFF3 file into an ``IntervalMetadata`` object
with upper bound of ``length``.

``DNA`` and ``Sequence`` GFF3 readers require ``seq_num`` of int as
parameter. It specifies which GFF3 record to read from a GFF3 file
with annotations of multiple sequences in it.

Writer-specific Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^
``skip_subregion`` is a boolean parameter used by all the GFF3 writers. It
specifies whether you would like to write each non-contiguous
sub-region for a feature annotation. For example, if there is
interval feature for a gene with two exons in an ``IntervalMetadata``
object, it will write one line into the GFF3 file when ``skip_subregion`` is
``True`` and will write 3 lines (one for the gene and one for each
exon, respectively) when ``skip_subregion`` is ``False``. Default is ``True``.

In addition, ``IntervalMetadata`` GFF3 writer needs a parameter of
``seq_id``. It specify the sequence ID (column 1 in GFF3 file) that
the annotation belong to.

Examples
--------

Let's create a file stream with following data in GFF3 format:

>>> gff_str = """
... ##gff-version\\t3.2.1
... ##sequence-region\\tctg123\\t1\\t1497228
... ctg123\\t.\\tgene\\t1000\\t9000\\t.\\t+\\t0\\tID=gene00001;Name=EDEN
... ctg123\\t.\\tTF_binding_site\\t1000\\t1012\\t.\\t+\\t.\\tParent=gene00001
... ctg123\\t.\\tmRNA\\t1050\\t9000\\t.\\t+\\t.\\tID=mRNA00001;Parent=gene00001
... """
>>> import io
>>> from skbio.metadata import IntervalMetadata
>>> from skbio.io import read
>>> gff = io.StringIO(gff_str)

We can read it into

>>> im = read(gff, format='gff3', into=IntervalMetadata,
...           seq_id='ctg123', length=10000)
>>> im   # doctest: +SKIP
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


References
----------
.. [1] https://github.com/The-Sequence-Ontology/\
Specifications/blob/master/gff3.md

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
def _gff3_to_generator(fh, lengths):
    '''Parse the GFF3 into the existing IntervalMetadata

    Parameters
    ----------
    fh : file
        file handler

    Yields
    ------
    tuple
        str of seq id, IntervalMetadata
    '''
    for (seq_id, lines), length in zip(_yield_record(fh), lengths):
        yield seq_id, _parse_record(lines, length)


@gff3.writer(None)
def _generator_to_gff3(obj, fh, skip_subregion=True):
    '''Write list of IntervalMetadata into file.

    Parameters
    ----------
    obj : Iterable of (seq_id, IntervalMetadata)
    fh : file handler
    '''
    for seq_id, obj_i in obj:
        _serialize_interval_metadata(obj_i, seq_id, fh, skip_subregion)


@gff3.reader(Sequence)
def _gff3_to_sequence(fh, seq_num=1):
    ''''''
    seq_id, lines = _get_nth_record(_yield_record(fh), seq_num)
    # you can't read directly from fh because in registry.py line 543
    # file.tell() will fail "telling position disabled by next() call".
    stream = StringIO(fh.read())
    seq = read(stream, format='fasta', into=Sequence, seq_num=seq_num)
    seq.interval_metadata = _parse_record(lines, len(seq))
    return seq


@gff3.writer(Sequence)
def _sequence_to_gff3(obj, fh, skip_subregion=True):
    _serialize_seq(obj, fh, skip_subregion)


@gff3.reader(DNA)
def _gff3_to_dna(fh, seq_num=1):
    ''''''
    seq_id, lines = _get_nth_record(_yield_record(fh), seq_num)
    stream = StringIO(fh.read())
    seq = read(stream, format='fasta', into=DNA, seq_num=seq_num)
    seq.interval_metadata = _parse_record(lines, len(seq))
    return seq


@gff3.writer(DNA)
def _dna_to_gff3(obj, fh, skip_subregion=True):
    _serialize_seq(obj, fh, skip_subregion)


@gff3.reader(IntervalMetadata)
def _gff3_to_interval_metadata(fh, seq_id, length):
    '''Read a GFF3 record into the specified interval metadata.

    Parameters
    ----------
    fh : file handler
    length : int
        seq length
    '''
    for sid, lines in _yield_record(fh):
        if sid == seq_id:
            return _parse_record(lines, length)
    # return an empty instead of None
    return IntervalMetadata(length)


@gff3.writer(IntervalMetadata)
def _interval_metadata_to_gff3(obj, fh, seq_id, skip_subregion=True):
    '''
    Parameters
    ----------
    obj : IntervalMetadata
    seq_id : str
        ID for column 1 in the GFF3 file.
    '''
    _serialize_interval_metadata(obj, seq_id, fh, skip_subregion=True)


def _yield_record(fh):
    '''Yield (seq_id, lines) that belong to the same sequence.'''
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


def _parse_record(lines, length):
    '''Parse the lines into a IntervalMetadata object.'''
    interval_metadata = IntervalMetadata(length)
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
        interval_metadata, seq_id, fh, skip_subregion=True):
    '''Serialize an IntervalMetadata to GFF3.

    Parameters
    ----------
    interval_metadata : IntervalMetadata
    skip_subregion : bool
        whether to skip outputting each sub region as a line in GFF3.
    '''
    # write file header
    fh.write('##gff-version 3\n')

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
        if len(bd) > 1 and skip_subregion is False:
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


def _serialize_seq(seq, fh, skip_subregion=True):
    '''Serialize a sequence to GFF3.'''
    _serialize_interval_metadata(
        seq.interval_metadata, seq.metadata['id'], fh, skip_subregion)
    print('##FASTA', file=fh)
    write(seq, into=fh, format='fasta')
