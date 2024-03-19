r"""QSeq format (:mod:`skbio.io.format.qseq`)
=========================================

.. currentmodule:: skbio.io.format.qseq

The QSeq format (`qseq`) is a record-based, plain text output format produced
by some DNA sequencers for storing biological sequence data, quality scores,
per-sequence filtering information, and run-specific metadata.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |generator of :mod:`skbio.sequence.Sequence` objects            |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
A QSeq file is composed of single-line records, delimited by tabs. There are
11 fields in a record:

- Machine name
- Run number
- Lane number (positive int)
- Tile number (positive int)
- X coordinate (integer)
- Y coordinate (integer)
- Index
- Read number (1-3)
- Sequence data (typically IUPAC characters)
- Quality scores (quality scores encoded as printable ASCII)
- Filter boolean (1 if sequence has passed CASAVA's filter, 0 otherwise)

For more details please refer to the CASAVA documentation [1]_.

.. note:: When a QSeq file is read into a scikit-bio object, the object's
   `metadata` attribute is automatically populated with data corresponding
   to the names above.

.. note:: `lowercase` functionality is supported when reading QSeq files.
   Refer to specific object constructor documentation for details.

.. note:: scikit-bio allows for the filter field to be ommitted, but it is not
   clear if this is part of the original format specification.

Format Parameters
-----------------
The following parameters are the same as in FASTQ format
(:mod:`skbio.io.format.fastq`):

- ``variant``: see ``variant`` parameter in FASTQ format
- ``phred_offset``: see ``phred_offset`` parameter in FASTQ format

The following additional parameters are the same as in FASTA format
(:mod:`skbio.io.format.fasta`):

- ``constructor``: see ``constructor`` parameter in FASTA format
- ``seq_num``: see ``seq_num`` parameter in FASTA format

Generators Only
^^^^^^^^^^^^^^^
- ``filter``: If `True`, excludes sequences that did not pass filtering
  (i.e., filter field is 0). Default is `True`.

Examples
--------
Suppose we have the following QSeq file::

    illumina	1	3	34	-30	30	0	1	ACG....ACGTAC	ruBBBBrBCEFGH	1
    illumina	1	3	34	30	-30	0	1	CGGGCATTGCA	CGGGCasdGCA	0
    illumina	1	3	35	-30	30	0	2	ACGTA.AATAAAC	geTaAafhwqAAf	1
    illumina	1	3	35	30	-30	0	3	CATTTAGGA.TGCA	tjflkAFnkKghvM	0

Let's define this file in-memory as a ``StringIO``, though this could be a real
file path, file handle, or anything that's supported by scikit-bio's I/O
registry in practice:

>>> from io import StringIO
>>> fs = '\n'.join([
...     'illumina\t1\t3\t34\t-30\t30\t0\t1\tACG....ACGTAC\truBBBBrBCEFGH\t1',
...     'illumina\t1\t3\t34\t30\t-30\t0\t1\tCGGGCATTGCA\tCGGGCasdGCA\t0',
...     'illumina\t1\t3\t35\t-30\t30\t0\t2\tACGTA.AATAAAC\tgeTaAafhwqAAf\t1',
...     'illumina\t1\t3\t35\t30\t-30\t0\t3\tCATTTAGGA.TGCA\ttjflkAFnkKghvM\t0'
... ])
>>> fh = StringIO(fs)

To iterate over the sequences using the generator reader, we run:

>>> import skbio.io
>>> for seq in skbio.io.read(fh, format='qseq', variant='illumina1.3'):
...     seq
...     print('')
Sequence
--------------------------------------
Metadata:
    'id': 'illumina_1:3:34:-30:30#0/1'
    'index': 0
    'lane_number': 3
    'machine_name': 'illumina'
    'read_number': 1
    'run_number': 1
    'tile_number': 34
    'x': -30
    'y': 30
Positional metadata:
    'quality': <dtype: uint8>
Stats:
    length: 13
--------------------------------------
0 ACG....ACG TAC
<BLANKLINE>
Sequence
--------------------------------------
Metadata:
    'id': 'illumina_1:3:35:-30:30#0/2'
    'index': 0
    'lane_number': 3
    'machine_name': 'illumina'
    'read_number': 2
    'run_number': 1
    'tile_number': 35
    'x': -30
    'y': 30
Positional metadata:
    'quality': <dtype: uint8>
Stats:
    length: 13
--------------------------------------
0 ACGTA.AATA AAC
<BLANKLINE>

Note that only two sequences were loaded because the QSeq reader filters out
sequences whose filter field is 0 (unless ``filter=False`` is supplied).

References
----------
.. [1] http://biowulf.nih.gov/apps/CASAVA_UG_15011196B.pdf


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.io import create_format, QSeqFormatError
from skbio.io.format._base import _decode_qual_to_phred, _get_nth_sequence
from skbio.sequence import Sequence, DNA, RNA, Protein

_default_phred_offset = None
_default_variant = None
_will_filter = True

qseq = create_format("qseq")


@qseq.sniffer()
def _qseq_sniffer(fh):
    empty = True
    try:
        for _, line in zip(range(10), fh):
            _record_parser(line)
            empty = False
        return not empty, {}
    except QSeqFormatError:
        return False, {}


@qseq.reader(None)
def _qseq_to_generator(
    fh,
    constructor=Sequence,
    filter=_will_filter,
    phred_offset=_default_phred_offset,
    variant=_default_variant,
    **kwargs,
):
    for line in fh:
        (
            machine_name,
            run,
            lane,
            tile,
            x,
            y,
            index,
            read,
            seq,
            raw_qual,
            filtered,
        ) = _record_parser(line)
        if not filter or not filtered:
            phred = _decode_qual_to_phred(raw_qual, variant, phred_offset)
            seq_id = "%s_%s:%s:%s:%s:%s#%s/%s" % (
                machine_name,
                run,
                lane,
                tile,
                x,
                y,
                index,
                read,
            )
            yield constructor(
                seq,
                metadata={
                    "id": seq_id,
                    "machine_name": machine_name,
                    "run_number": int(run),
                    "lane_number": int(lane),
                    "tile_number": int(tile),
                    "x": int(x),
                    "y": int(y),
                    "index": int(index),
                    "read_number": int(read),
                },
                positional_metadata={"quality": phred},
                **kwargs,
            )


@qseq.reader(Sequence)
def _qseq_to_sequence(
    fh,
    seq_num=1,
    phred_offset=_default_phred_offset,
    variant=_default_variant,
    **kwargs,
):
    return _get_nth_sequence(
        _qseq_to_generator(
            fh,
            filter=False,
            phred_offset=phred_offset,
            variant=variant,
            constructor=Sequence,
            **kwargs,
        ),
        seq_num,
    )


@qseq.reader(DNA)
def _qseq_to_dna(
    fh,
    seq_num=1,
    phred_offset=_default_phred_offset,
    variant=_default_variant,
    **kwargs,
):
    return _get_nth_sequence(
        _qseq_to_generator(
            fh,
            filter=False,
            phred_offset=phred_offset,
            variant=variant,
            constructor=DNA,
            **kwargs,
        ),
        seq_num,
    )


@qseq.reader(RNA)
def _qseq_to_rna(
    fh,
    seq_num=1,
    phred_offset=_default_phred_offset,
    variant=_default_variant,
    **kwargs,
):
    return _get_nth_sequence(
        _qseq_to_generator(
            fh,
            filter=False,
            phred_offset=phred_offset,
            variant=variant,
            constructor=RNA,
            **kwargs,
        ),
        seq_num,
    )


@qseq.reader(Protein)
def _qseq_to_protein(
    fh,
    seq_num=1,
    phred_offset=_default_phred_offset,
    variant=_default_variant,
    **kwargs,
):
    return _get_nth_sequence(
        _qseq_to_generator(
            fh,
            filter=False,
            phred_offset=phred_offset,
            variant=variant,
            constructor=Protein,
            **kwargs,
        ),
        seq_num,
    )


def _record_parser(line):
    fields = line.rstrip("\n")
    if fields:
        fields = fields.split("\t")
    else:
        raise QSeqFormatError("Found blank line.")
    f_len = len(fields)
    if not (10 <= f_len <= 11):
        raise QSeqFormatError("Expected 10 or 11 fields, found %d." % f_len)
    # If the filter field was ommitted, assume that it passed filtering:
    if f_len == 10:
        fields.append("1")

    (machine, run, lane, tile, x, y, index, read, seq, raw_qaul, filter) = fields

    _test_fields([("filter", filter)], lambda x: x in "01", "0 or 1")

    _test_fields([("read", read)], lambda x: x in "123", "in the range [1, 3]")

    _test_fields([("x", x), ("y", y)], lambda x: int(x) is not None, "an integer")

    _test_fields(
        [("lane", lane), ("tile", tile)], lambda x: int(x) >= 0, "a positive integer"
    )

    return (machine, run, lane, tile, x, y, index, read, seq, raw_qaul, filter == "0")


def _test_fields(iterkv, test, efrag):
    try:
        for k, v in iterkv:
            if not test(v):
                raise ValueError()
    except ValueError:
        raise QSeqFormatError("Field %r is not %s." % (k, efrag))
