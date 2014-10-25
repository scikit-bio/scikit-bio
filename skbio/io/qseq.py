"""
Newick format (:mod:`skbio.io.qseq`)
======================================

.. currentmodule:: skbio.io.qseq

Format Support
--------------
**Has Sniffer: Yes**

+----------+----------+------------------------------------------------------+
|**Reader**|**Writer**|                   **Object Class**                   |
+----------+----------+------------------------------------------------------+
|Yes       |No        |generator of :mod:`skbio.sequence.BiologicalSequence` |
|          |          |objects                                               |
+----------+----------+------------------------------------------------------+
|Yes       |No        |:mod:`skbio.alignment.SequenceCollection`             |
+----------+----------+------------------------------------------------------+
|Yes       |No        |:mod:`skbio.sequence.BiologicalSequence`              |
+----------+----------+------------------------------------------------------+
|Yes       |No        |:mod:`skbio.sequence.NucleotideSequence`              |
+----------+----------+------------------------------------------------------+
|Yes       |No        |:mod:`skbio.sequence.DNASequence`                     |
+----------+----------+------------------------------------------------------+
|Yes       |No        |:mod:`skbio.sequence.RNASequence`                     |
+----------+----------+------------------------------------------------------+
|Yes       |No        |:mod:`skbio.sequence.ProteinSequence`                 |
+----------+----------+------------------------------------------------------+

Format Specification
--------------------

Examples
--------


References
----------
.. [1] http://evolution.genetics.washington.edu/phylip/newick_doc.html
.. [2] http://evolution.genetics.washington.edu/phylip/newicktree.html

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
    if f_len == 10:
        fields.append('1')

    (machine, run, lane, tile, x, y, index, read, seq, raw_qaul,
     filter) = tuple(fields)

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
