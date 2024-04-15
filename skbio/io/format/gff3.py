"""GFF3 format (:mod:`skbio.io.format.gff3`)
=========================================

.. currentmodule:: skbio.io.format.gff3

GFF3 (Generic Feature Format version 3) is a standard file format for
describing features for biological sequences. It contains lines of
text, each consisting of 9 tab-delimited columns [1]_.

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
``IntervalMetadata`` GFF3 reader requires 1 parameter: ``seq_id``.
It reads the annotation with the specified
sequence ID from the GFF3 file into an ``IntervalMetadata`` object.

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

>>> from skbio import Sequence, DNA
>>> gff_str = '''
... ##gff-version 3
... seq_1\\t.\\tgene\\t10\\t90\\t.\\t+\\t0\\tID=gen1
... seq_1\\t.\\texon\\t10\\t30\\t.\\t+\\t.\\tParent=gen1
... seq_1\\t.\\texon\\t50\\t90\\t.\\t+\\t.\\tParent=gen1
... seq_2\\t.\\tgene\\t80\\t96\\t.\\t-\\t.\\tID=gen2
... ##FASTA
... >seq_1
... ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
... ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
... >seq_2
... ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
... ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
... '''
>>> import io
>>> from skbio.metadata import IntervalMetadata
>>> from skbio.io import read
>>> gff = io.StringIO(gff_str)

We can read it into ``IntervalMetadata``. Each line will be read into
an interval feature in ``IntervalMetadata`` object:

>>> im = read(gff, format='gff3', into=IntervalMetadata,
...           seq_id='seq_1')
>>> im  # doctest: +SKIP
3 interval features
-------------------
Interval(interval_metadata=<4604421736>, bounds=[(9, 90)], \
fuzzy=[(False, False)], metadata={'type': 'gene', \
'phase': 0, 'strand': '+', 'source': '.', 'score': '.', 'ID': 'gen1'})
Interval(interval_metadata=<4604421736>, bounds=[(9, 30)], \
fuzzy=[(False, False)], metadata={'strand': '+', 'source': '.', \
'type': 'exon', 'Parent': 'gen1', 'score': '.'})
Interval(interval_metadata=<4604421736>, bounds=[(49, 90)], \
fuzzy=[(False, False)], metadata={'strand': '+', 'source': '.', \
'type': 'exon', 'Parent': 'gen1', 'score': '.'})

We can write the ``IntervalMetadata`` object back to GFF3 file:

>>> with io.StringIO() as fh:    # doctest: +NORMALIZE_WHITESPACE
...     print(im.write(fh, format='gff3', seq_id='seq_1').getvalue())
##gff-version 3
seq_1	.	gene	10	90	.	+	0	ID=gen1
seq_1	.	exon	10	30	.	+	.	Parent=gen1
seq_1	.	exon	50	90	.	+	.	Parent=gen1
<BLANKLINE>

If the GFF3 file does not have the sequence ID, it will return an empty object:

>>> gff = io.StringIO(gff_str)
>>> im = read(gff, format='gff3', into=IntervalMetadata,
...           seq_id='foo')
>>> im
0 interval features
-------------------

We can also read the GFF3 file into a generator:

>>> gff = io.StringIO(gff_str)
>>> gen = read(gff, format='gff3')
>>> for im in gen:   # doctest: +SKIP
...     print(im[0])   # the seq id
...     print(im[1])   # the interval metadata on this seq
seq_1
3 interval features
-------------------
Interval(interval_metadata=<4603377592>, bounds=[(9, 90)], \
fuzzy=[(False, False)], metadata={'type': 'gene', 'ID': 'gen1', \
'source': '.', 'score': '.', 'strand': '+', 'phase': 0})
Interval(interval_metadata=<4603377592>, bounds=[(9, 30)], \
fuzzy=[(False, False)], metadata={'strand': '+', 'type': 'exon', \
'Parent': 'gen1', 'source': '.', 'score': '.'})
Interval(interval_metadata=<4603377592>, bounds=[(49, 90)], \
fuzzy=[(False, False)], metadata={'strand': '+', 'type': 'exon', \
'Parent': 'gen1', 'source': '.', 'score': '.'})
seq_2
1 interval feature
------------------
Interval(interval_metadata=<4603378712>, bounds=[(79, 96)], \
fuzzy=[(False, False)], metadata={'strand': '-', 'type': 'gene', \
'ID': 'gen2', 'source': '.', 'score': '.'})

For the GFF3 file with sequences, we can read it into ``Sequence`` or ``DNA``:

>>> gff = io.StringIO(gff_str)
>>> seq = read(gff, format='gff3', into=Sequence, seq_num=1)
>>> seq
Sequence
--------------------------------------------------------------------
Metadata:
    'description': ''
    'id': 'seq_1'
Interval metadata:
    3 interval features
Stats:
    length: 100
--------------------------------------------------------------------
0  ATGCATGCAT GCATGCATGC ATGCATGCAT GCATGCATGC ATGCATGCAT GCATGCATGC
60 ATGCATGCAT GCATGCATGC ATGCATGCAT GCATGCATGC

>>> gff = io.StringIO(gff_str)
>>> seq = read(gff, format='gff3', into=DNA, seq_num=2)
>>> seq
DNA
--------------------------------------------------------------------
Metadata:
    'description': ''
    'id': 'seq_2'
Interval metadata:
    1 interval feature
Stats:
    length: 120
    has gaps: False
    has degenerates: False
    has definites: True
    GC-content: 50.00%
--------------------------------------------------------------------
0  ATGCATGCAT GCATGCATGC ATGCATGCAT GCATGCATGC ATGCATGCAT GCATGCATGC
60 ATGCATGCAT GCATGCATGC ATGCATGCAT GCATGCATGC ATGCATGCAT GCATGCATGC


References
----------
.. [1] https://github.com/The-Sequence-Ontology/\
Specifications/blob/master/gff3.md


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import re
from collections.abc import Iterable

from skbio.sequence import DNA, Sequence
from skbio.io import create_format, GFF3FormatError
from skbio.metadata import IntervalMetadata
from skbio.io.format._base import _line_generator, _too_many_blanks, _get_nth_sequence
from skbio.io.format.fasta import _fasta_to_generator
from skbio.io.format._sequence_feature_vocabulary import (
    _vocabulary_change,
    _vocabulary_skip,
)
from skbio.io import write


gff3 = create_format("gff3")


@gff3.sniffer()
def _gff3_sniffer(fh):
    # check the 1st real line is a valid ID line
    if _too_many_blanks(fh, 5):
        return False, {}

    try:
        line = next(_line_generator(fh, skip_blanks=True, strip=False))
    except StopIteration:
        return False, {}

    if re.match(r"##gff-version\s+3", line):
        return True, {}
    else:
        return False, {}


@gff3.reader(None)
def _gff3_to_generator(fh):
    """Parse the GFF3 into the existing IntervalMetadata.

    Parameters
    ----------
    fh : file
        file handler

    Yields
    ------
    tuple
        str of seq id, IntervalMetadata

    """
    id_lengths = {}
    for data_type, sid, data in _yield_record(fh):
        if data_type == "length":
            # get length from sequence-region pragma.
            # the pragma lines are always before the real annotation lines.
            id_lengths[sid] = data
        elif data_type == "data":
            length = id_lengths.get(sid)
            yield sid, _parse_record(data, length)


@gff3.writer(None)
def _generator_to_gff3(obj, fh, skip_subregion=True):
    """Write a list of IntervalMetadata into a file.

    Parameters
    ----------
    obj : Iterable of (seq_id, IntervalMetadata)
        List of IntervalMetadata to write.
    fh : file handler
        File to write into.
    skip_subregion : bool
        Write a line for each sub-regions of an ``Interval`` if it is ``False``.

    """
    # write file header
    fh.write("##gff-version 3\n")
    for seq_id, obj_i in obj:
        _serialize_interval_metadata(obj_i, seq_id, fh, skip_subregion)


@gff3.reader(Sequence)
def _gff3_to_sequence(fh, seq_num=1):
    return _construct_seq(fh, Sequence, seq_num)


@gff3.writer(Sequence)
def _sequence_to_gff3(obj, fh, skip_subregion=True):
    # write file header
    fh.write("##gff-version 3\n")
    _serialize_seq(obj, fh, skip_subregion)


@gff3.reader(DNA)
def _gff3_to_dna(fh, seq_num=1):
    return _construct_seq(fh, DNA, seq_num)


@gff3.writer(DNA)
def _dna_to_gff3(obj, fh, skip_subregion=True):
    # write file header
    fh.write("##gff-version 3\n")
    _serialize_seq(obj, fh, skip_subregion)


@gff3.reader(IntervalMetadata)
def _gff3_to_interval_metadata(fh, seq_id):
    """Read a GFF3 record into the specified interval metadata.

    Parameters
    ----------
    fh : file handler
        GFF3 file to read.
    seq_id : str
        Sequence ID which the interval metadata is associated with.

    """
    length = None
    for data_type, sid, data in _yield_record(fh):
        if seq_id == sid:
            if data_type == "length":
                # get length from sequence-region pragma
                length = data
            elif data_type == "data":
                return _parse_record(data, length)
            else:
                raise GFF3FormatError(
                    "Unknown section in the input GFF3 file: "
                    "%r %r %r" % (data_type, sid, data)
                )
    # return an empty instead of None
    return IntervalMetadata(None)


@gff3.writer(IntervalMetadata)
def _interval_metadata_to_gff3(obj, fh, seq_id, skip_subregion=True):
    """Output ``IntervalMetadata`` object to GFF3 file.

    Parameters
    ----------
    obj : IntervalMetadata
        An IntervalMetadata object.
    fh : file object
        File object opened for writing.
    seq_id : str
        ID for column 1 in the GFF3 file.
    skip_subregion : bool
        Write a line for each sub-regions of an ``Interval`` if it is ``False``.

    """
    # write file header
    fh.write("##gff-version 3\n")
    _serialize_interval_metadata(obj, seq_id, fh, skip_subregion=True)


def _construct_seq(fh, constructor=DNA, seq_num=1):
    lines = []
    for i, (data_type, seq_id, L) in enumerate(_yield_record(fh), 1):
        if data_type == "data" and seq_num == i:
            lines = L
    seq = _get_nth_sequence(
        _fasta_to_generator(fh, constructor=constructor), seq_num=seq_num
    )
    seq.interval_metadata = _parse_record(lines, len(seq))
    return seq


def _yield_record(fh):
    """Yield (seq_id, lines) that belong to the same sequence."""
    lines = []
    current = False
    for line in _line_generator(fh, skip_blanks=True, strip=True):
        if line.startswith("##sequence-region"):
            _, seq_id, start, end = line.split()
            length = int(end) - int(start) + 1
            yield "length", seq_id, length
        if line.startswith("##FASTA"):
            # stop once reaching to sequence section
            break
        if not line.startswith("#"):
            try:
                seq_id, _ = line.split("\t", 1)
            except ValueError:
                raise GFF3FormatError("Wrong GFF3 format at line: %s" % line)
            if current == seq_id:
                lines.append(line)
            else:
                if current is not False:
                    yield "data", current, lines
                lines = [line]
                current = seq_id
    if current is False:
        # if the input file object is empty, it should return
        # an empty generator
        return
        yield
    else:
        yield "data", current, lines


def _parse_record(lines, length):
    """Parse the lines into a IntervalMetadata object."""
    interval_metadata = IntervalMetadata(length)
    for line in lines:
        columns = line.split("\t")
        # there should be 9 columns
        if len(columns) != 9:
            raise GFF3FormatError('do not have 9 columns in this line: "%s"' % line)
        # the 1st column is seq ID for every feature. don't store
        # this repetitive information
        metadata = {
            "source": columns[1],
            "type": columns[2],
            "score": columns[5],
            "strand": columns[6],
        }
        phase = columns[7]
        # phase value can only be int or '.'
        try:
            metadata["phase"] = int(phase)
        except ValueError:
            if phase != ".":
                raise GFF3FormatError(
                    "unknown value for phase column: {!r}".format(phase)
                )
        metadata.update(_parse_attr(columns[8]))

        start, end = columns[3:5]

        bounds = [(int(start) - 1, int(end))]

        interval_metadata.add(bounds, metadata=metadata)

    return interval_metadata


def _parse_attr(s):
    """Parse attribute column."""
    voca_change = _vocabulary_change("gff3")
    md = {}
    # in case the line ending with ';', strip it.
    s = s.rstrip(";")
    for attr in s.split(";"):
        k, v = attr.split("=")
        if k in voca_change:
            k = voca_change[k]
        md[k] = v
    return md


def _serialize_interval_metadata(interval_metadata, seq_id, fh, skip_subregion=True):
    """Serialize an IntervalMetadata to GFF3.

    Parameters
    ----------
    interval_metadata : IntervalMetadata
        An IntervalMetadata object.
    seq_id : str
        Seq id for the current annotation. It will be used as the 1st column
        in the GFF3.
    fh : file handler
        The file object to output.
    skip_subregion : bool
        Whether to skip outputting each sub region as a line in GFF3.

    """
    column_keys = ["source", "type", "score", "strand", "phase"]
    voca_change = _vocabulary_change("gff3", False)
    voca_skip = _vocabulary_skip("gff3")
    voca_skip.extend(column_keys)

    # these characters have reserved meanings in column 9 and must be
    # escaped when used in other contexts
    escape = str.maketrans({";": "%3B", "=": "%3D", "&": "%26", ",": "%2C"})

    for interval in interval_metadata._intervals:
        md = interval.metadata
        bd = interval.bounds
        start = str(bd[0][0] + 1)
        end = str(bd[-1][-1])

        source, feat_type, score, strand, phase = [
            str(md.get(i, ".")) for i in column_keys
        ]
        columns = [seq_id, source, feat_type, start, end, score, strand, phase]

        # serialize the attributes in column 9
        attr = []
        # use sort to make the output order deterministic
        for k in sorted(md):
            if k in voca_skip:
                # skip the metadata that doesn't go to attribute column
                continue
            v = md[k]
            if k in voca_change:
                k = voca_change[k]
            if isinstance(v, Iterable) and not isinstance(v, str):
                # if there are multiple values for this attribute,
                # convert them to str and concat them with ","
                v = ",".join(str(i).translate(escape) for i in v)
            else:
                v = v.translate(escape)
            attr.append("%s=%s" % (k.translate(escape), v))

        columns.append(";".join(attr))

        fh.write("\t".join(columns))
        fh.write("\n")

        # if there are multiple regions for this feature,
        # output each region as a standalone line in GFF3.
        if len(bd) > 1 and skip_subregion is False:
            for start, end in bd:
                # if this is a gene, then each sub region should be an exon
                if columns[2] == "gene":
                    columns[2] = "exon"
                columns[3] = str(start + 1)
                columns[4] = str(end)
                try:
                    parent = md["ID"]
                except KeyError:
                    raise GFF3FormatError(
                        "You need provide ID info for "
                        "the parent interval feature: %r" % interval
                    )
                columns[8] = "Parent=%s" % parent
                fh.write("\t".join(columns))
                fh.write("\n")


def _serialize_seq(seq, fh, skip_subregion=True):
    """Serialize a sequence to GFF3."""
    _serialize_interval_metadata(
        seq.interval_metadata, seq.metadata["id"], fh, skip_subregion
    )
    fh.write("##FASTA\n")
    write(seq, into=fh, format="fasta")
