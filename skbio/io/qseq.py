r"""
QSeq format (:mod:`skbio.io.qseq`)
==================================

.. currentmodule:: skbio.io.qseq

The QSeq format (`qseq`) is a record-based, plain text output format produced
by some DNA sequencers for storing biological sequence data, quality scores,
per-sequence filtering information, and run-specific metadata.

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |No    |generator of :mod:`skbio.sequence.BiologicalSequence` objects  |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.alignment.SequenceCollection`                      |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.BiologicalSequence`                       |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.NucleotideSequence`                       |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.DNASequence`                              |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.RNASequence`                              |
+------+------+---------------------------------------------------------------+
|Yes   |No    |:mod:`skbio.sequence.ProteinSequence`                          |
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

.. note:: scikit-bio allows for the filter field to be ommitted, but it is not
   clear if this is part of the original format specification.

Format Parameters
-----------------
The following parameters are the same as in FASTQ format
(:mod:`skbio.io.fastq`):

- ``variant``: see ``variant`` parameter in FASTQ format
- ``phred_offset``: see ``phred_offset`` parameter in FASTQ format

The following additional parameters are the same as in FASTA format
(:mod:`skbio.io.fasta`):

- ``constructor``: see ``constructor`` parameter in FASTA format
- ``seq_num``: see ``seq_num`` parameter in FASTA format

SequenceCollection and Generators Only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

>>> from StringIO import StringIO
>>> fs = '\n'.join([
...     'illumina\t1\t3\t34\t-30\t30\t0\t1\tACG....ACGTAC\truBBBBrBCEFGH\t1',
...     'illumina\t1\t3\t34\t30\t-30\t0\t1\tCGGGCATTGCA\tCGGGCasdGCA\t0',
...     'illumina\t1\t3\t35\t-30\t30\t0\t2\tACGTA.AATAAAC\tgeTaAafhwqAAf\t1',
...     'illumina\t1\t3\t35\t30\t-30\t0\t3\tCATTTAGGA.TGCA\ttjflkAFnkKghvM\t0'
... ])
>>> fh = StringIO(fs)

To load the sequences into a ``SequenceCollection``, we run:

>>> from skbio import SequenceCollection
>>> sc = SequenceCollection.read(fh, variant='illumina1.3')
>>> sc
<SequenceCollection: n=2; mean +/- std length=13.00 +/- 0.00>

Note that only two sequences were loaded because the QSeq reader filters out
sequences whose filter field is 0 (unless ``filter=False`` is supplied).

References
----------
.. [1] http://biowulf.nih.gov/apps/CASAVA_UG_15011196B.pdf

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from future.builtins import zip, range

from skbio.io import register_reader, register_sniffer, QSeqFormatError
from skbio.io._base import _decode_qual_to_phred, _get_nth_sequence
from skbio.alignment import SequenceCollection
from skbio.sequence import (BiologicalSequence, NucleotideSequence,
                            DNASequence, RNASequence, ProteinSequence)

_default_phred_offset = None
_default_variant = None
_will_filter = True


@register_sniffer('qseq')
def _qseq_sniffer(fh):
    empty = True
    try:
        for _, line in zip(range(10), fh):
            _record_parser(line)
            empty = False
        return not empty, {}
    except QSeqFormatError:
        return False, {}


@register_reader('qseq')
def _qseq_to_generator(fh, constructor=BiologicalSequence, filter=_will_filter,
                       phred_offset=_default_phred_offset,
                       variant=_default_variant):
    for line in fh:
        (machine_name, run, lane, tile, x, y, index, read, seq, raw_qual,
         filtered) = _record_parser(line)
        if not filter or not filtered:
            phred = _decode_qual_to_phred(raw_qual, variant, phred_offset)
            seq_id = '%s_%s:%s:%s:%s:%s#%s/%s' % (
                machine_name, run, lane, tile, x, y, index, read)
            yield constructor(seq, quality=phred, id=seq_id)


@register_reader('qseq', SequenceCollection)
def _qseq_to_sequence_collection(fh, constructor=BiologicalSequence,
                                 filter=_will_filter,
                                 phred_offset=_default_phred_offset,
                                 variant=_default_variant):
    return SequenceCollection(list(_qseq_to_generator(
        fh, constructor=constructor, filter=filter, phred_offset=phred_offset,
        variant=variant)))


@register_reader('qseq', BiologicalSequence)
def _qseq_to_biological_sequence(fh, seq_num=1,
                                 phred_offset=_default_phred_offset,
                                 variant=_default_variant):
    return _get_nth_sequence(_qseq_to_generator(fh, filter=False,
                             phred_offset=phred_offset, variant=variant,
                             constructor=BiologicalSequence), seq_num)


@register_reader('qseq', NucleotideSequence)
def _qseq_to_nucleotide_sequence(fh, seq_num=1,
                                 phred_offset=_default_phred_offset,
                                 variant=_default_variant):
    return _get_nth_sequence(_qseq_to_generator(fh, filter=False,
                             phred_offset=phred_offset, variant=variant,
                             constructor=NucleotideSequence), seq_num)


@register_reader('qseq', DNASequence)
def _qseq_to_dna_sequence(fh, seq_num=1,
                          phred_offset=_default_phred_offset,
                          variant=_default_variant):
    return _get_nth_sequence(_qseq_to_generator(fh, filter=False,
                             phred_offset=phred_offset, variant=variant,
                             constructor=DNASequence), seq_num)


@register_reader('qseq', RNASequence)
def _qseq_to_rna_sequence(fh, seq_num=1,
                          phred_offset=_default_phred_offset,
                          variant=_default_variant):
    return _get_nth_sequence(_qseq_to_generator(fh, filter=False,
                             phred_offset=phred_offset, variant=variant,
                             constructor=RNASequence), seq_num)


@register_reader('qseq', ProteinSequence)
def _qseq_to_protein_sequence(fh, seq_num=1,
                              phred_offset=_default_phred_offset,
                              variant=_default_variant):
    return _get_nth_sequence(_qseq_to_generator(fh, filter=False,
                             phred_offset=phred_offset, variant=variant,
                             constructor=ProteinSequence), seq_num)


def _record_parser(line):
    fields = line.rstrip('\n')
    if fields:
        fields = fields.split('\t')
    else:
        raise QSeqFormatError('Found blank line.')
    f_len = len(fields)
    if not (10 <= f_len <= 11):
        raise QSeqFormatError('Expected 10 or 11 fields, found %d.' % f_len)
    # If the filter field was ommitted, assume that it passed filtering:
    if f_len == 10:
        fields.append('1')

    (machine, run, lane, tile, x, y, index, read, seq, raw_qaul,
     filter) = fields

    _test_fields([('filter', filter)], lambda x: x in '01',
                 "0 or 1")

    _test_fields([('read', read)], lambda x: x in '123',
                 "in the range [1, 3]")

    _test_fields([('x', x), ('y', y)], lambda x: int(x) is not None,
                 "an integer")

    _test_fields([('lane', lane), ('tile', tile)], lambda x: int(x) >= 0,
                 "a positive integer")

    return (machine, run, lane, tile, x, y, index, read, seq, raw_qaul,
            filter == '0')


def _test_fields(iterkv, test, efrag):
    try:
        for k, v in iterkv:
            if not test(v):
                raise ValueError()
    except ValueError:
        raise QSeqFormatError('Field %r is not %s.' % (k, efrag))
