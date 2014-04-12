#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from itertools import chain, izip
from collections import namedtuple

from ..parse.sequences import parse_fasta, parse_fastq, parse_qual
from .workflow import Workflow, not_none, method, requires


class SequenceRecord(namedtuple("SequenceRecord", ("SequenceID", "Sequence",
                                                   "QualID", "Qual"))):
    __slots__ = ()

    def __new__(cls, SequenceID, Sequence, QualID=None, Qual=None):
        return super(SequenceRecord, cls).__new__(cls, SequenceID, Sequence,
                                                  QualID, Qual)


def _has_qual(item):
    """Return True if it appears that there is qual data"""
    return (item['QualID'] is not None) and (item['Qual'] is not None)


class SequenceIterator(Workflow):
    """Provide a standard API for interacting with sequence files

    Provide a common interface for iterating over sequence files, including
    support for quality scores and transforms.

    A transform method is a function that takes the state dict and modifies it
    in place. For instance, to reverse sequences, you could pass in the
    following function:

    def reverse(st):
        st['Sequence']= st['Sequence'][::-1]
        st['Qual'] = st['Qual'][::-1] if st['Qual'] is not None else None

    as transform. The primary intention is to support reverse complementing
    of sequences.

    All subclasses of this object are expected to yield a dict that contains
    the following keys and types

        SequenceID : str, the sequence identifier
        Sequence   : str, the sequence itself
        QualID     : str or None, the quality ID (for completeness)
        Qual       : np.array or None, the quality scores

    The yielded object is preallocated a single time, as such (assuming the
    first two sequences are in fact not identical)::

        gen = instance_of_sequence_iterator()

        item1 = gen.next()
        seq1 = item1['Sequence']

        item2 = gen.next()
        seq2 = item2['Sequence']

        id(item1) == id(item2) # True
        item1 == item2         # True

        id(seq1) == id(seq2)   # False
        seq1 == seq2           # False
    """
    def __init__(self, seq, qual=None, transform=None, valid_id=True,
                 valid_length=True):
        self.seq = seq
        self.qual = qual

        self._transform = transform

        state = {'SequenceID': None,
                 'Sequence': None,
                 'QualID': None,
                 'Qual': None}

        options = {'transform': self._transform,
                   'valid_id': valid_id,
                   'valid_length': valid_length}

        if self.seq is None:
            raise ValueError("SequenceIterator requires sequences!")

        super(SequenceIterator, self).__init__(state, options=options)

    def _gen(self):
        """Yield a populated record"""
        raise NotImplementedError("Must be implemented by subclass")

    def __call__(self):
        return super(SequenceIterator, self).__call__(self._gen())

    def __iter__(self):
        return self()

    def initialize_state(self, item):
        """Reset the buffer"""
        self.state['SequenceID'] = item.SequenceID
        self.state['Sequence'] = item.Sequence
        self.state['QualID'] = item.QualID
        self.state['Qual'] = item.Qual

    @method(priority=100)
    @requires(option='valid_id', values=True, state=_has_qual)
    def validate_ids(self):
        self.failed = self.state['SequenceID'] != self.state['QualID']

    @method(priority=90)
    @requires(option='valid_length', values=True, state=_has_qual)
    def valid_lengths(self):
        self.failed = len(self.state['Sequence']) != len(self.state['Qual'])

    @method(priority=80)
    @requires(option='transform', values=not_none)
    def transform(self):
        """Transform state if necessary"""
        self._transform(self.state)


class FastaIterator(SequenceIterator):
    """Populate and yield records based on fasta sequence"""
    def _gen(self):
        """Construct internal iterators"""
        # construct fasta generators
        fasta_gens = chain(*[parse_fasta(f) for f in self.seq])

        # construct qual generators if necessary
        if self.qual is not None:
            qual_gens = chain(*[parse_qual(f) for f in self.qual])
        else:
            qual_gens = None

        # determine which specific generator to return
        if qual_gens is None:
            gen = self._fasta_gen(fasta_gens)
        else:
            gen = self._fasta_qual_gen(fasta_gens, qual_gens)

        return gen

    def _fasta_gen(self, fasta_gens):
        """Yield fasta data"""
        _iter = fasta_gens
        for (seq_id, seq) in _iter:
            yield SequenceRecord(seq_id, seq)

    def _fasta_qual_gen(self, fasta_gen, qual_gen):
        """Yield fasta and qual together"""
        _iter = izip(fasta_gen, qual_gen)
        for (seq_id, seq), (qual_id, qual) in _iter:
            yield SequenceRecord(seq_id, seq, qual_id, qual)


class FastqIterator(SequenceIterator):
    """Populate and yield records based on fastq sequence

    Note: thq 'qual' keyword argument is ignored by this object.
    """
    def _gen(self):
        """Construct internal iterators"""
        # construct fastq parsers
        fastq_gens = chain(*[parse_fastq(f, False) for f in self.seq])

        # peek the first record to determine phred offset
        first_item = next(fastq_gens)
        seqid, seq, qual = first_item
        fastq_gens = chain([first_item], fastq_gens)

        gen = self._fastq_gen(fastq_gens)

        return gen

    def _fastq_gen(self, fastq_gens):
        """Yield fastq data"""
        for (seq_id, seq, qual) in fastq_gens:
            yield SequenceRecord(seq_id, seq, seq_id, qual)
