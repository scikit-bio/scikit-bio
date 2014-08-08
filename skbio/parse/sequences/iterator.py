# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import chain

from future.builtins import zip

from skbio.core.workflow import Workflow, not_none, method, requires
from .fasta import parse_fasta, parse_qual
from .fastq import parse_fastq
from .qseq import parse_qseq


def _has_qual(item):
    """Return True if it appears that there is qual data"""
    return (item['QualID'] is not None) and (item['Qual'] is not None)


class SequenceIterator(Workflow):
    """Provide a standard API for interacting with sequence data

    Provide a common interface for iterating over sequence data, including
    support for quality scores and transforms.

    A transform method is a function that takes the state dict and modifies it
    in place. For instance, to reverse sequences, you could pass in the
    following function:

    def reverse(st):
        st['Sequence']= st['Sequence'][::-1]
        st['Qual'] = st['Qual'][::-1] if st['Qual'] is not None else None

    as ``transform``. The primary intention is to support reverse complementing
    of sequences.

    All subclasses of this object are expected to update the following in
    ``state``:

        SequenceID : str, the sequence identifier
        Sequence   : str, the sequence itself
        QualID     : str or None, the quality ID (for completeness)
        Qual       : np.array or None, the quality scores

    ``state`` is preallocated a single time to avoid repetitive allocations.
    What this means is that the object being yielded is updated in place. If
    an individual record needs to be tracked over time, then it is recommended
    that copies of the yielded data are made.

    *WARNING*: The yielded obj is not safe for use with Python 2.7's builtin
    `zip` method as the state is updated in place.

    Parameters
    ----------

    seq : list of open file-like objects
    qual : list of open file-like objects or None
    transform : function or None
        If provided, this function will be passed ``state``
    valid_id : bool
        If true, verify sequence and qual IDs are identical (if relevant)
    valid_length : bool
        If true, verify the length of the sequence and qual are the same
        (if relevant)

    Attributes
    ----------

    seq
    qual
    state
    options

    """
    def __init__(self, seq, qual=None, transform=None, valid_id=True,
                 valid_length=True, **kwargs):
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
        """Do nothing here as the subclassed iterators update state directly"""
        pass

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
    """Populate state based on fasta sequence and qual (if provided)"""
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
            self.state['SequenceID'] = seq_id
            self.state['Sequence'] = seq

            # as we're updating state in place and effectively circumventing
            # Workflow.initialize_state, we do not need to yield anything
            yield None

    def _fasta_qual_gen(self, fasta_gen, qual_gen):
        """Yield fasta and qual together"""
        _iter = zip(fasta_gen, qual_gen)
        for (seq_id, seq), (qual_id, qual) in _iter:
            self.state['SequenceID'] = seq_id
            self.state['Sequence'] = seq
            self.state['QualID'] = qual_id
            self.state['Qual'] = qual

            # as we're updating state in place and effectively circumventing
            # Workflow.initialize_state, we do not need to yield anything
            yield None


class FastqIterator(SequenceIterator):
    """Populate state based on fastq sequence

    Note: thq 'qual' keyword argument is ignored by this object.
    """
    def __init__(self, *args, **kwargs):
        if 'phred_offset' in kwargs:
            self._fpo = kwargs.pop('phred_offset')
        else:
            # force to an offset of 33
            self._fpo = 33

        super(FastqIterator, self).__init__(*args, **kwargs)

    def _gen(self):
        """Construct internal iterators"""
        fastq_gens = chain(*[parse_fastq(f, phred_offset=self._fpo)
                             for f in self.seq])
        return self._fastq_gen(fastq_gens)

    def _fastq_gen(self, fastq_gens):
        """Yield fastq data"""
        for (seq_id, seq, qual) in fastq_gens:
            self.state['SequenceID'] = seq_id
            self.state['Sequence'] = seq
            self.state['QualID'] = seq_id
            self.state['Qual'] = qual

            # as we're updating state in place and effectively circumventing
            # Workflow.initialize_state, we do not need to yield anything
            yield None


class QseqIterator(SequenceIterator):
    """Populate state based on qseq sequence."""
    def __init__(self, *args, **kwargs):
        if 'phred_offset' in kwargs:
            self._fpo = kwargs.pop('phred_offset')
        else:
            # force to an offset of 33
            self._fpo = 33

        super(QseqIterator, self).__init__(*args, **kwargs)

    def _gen(self):
        """Construct internal iterators"""
        # construct qseq generators
        qseq_gens = chain(*[parse_qseq(f,  phred_offset=self._fpo)
                            for f in self.seq])

        gen = self._qseq_gen(qseq_gens)

        return gen

    def _qseq_gen(self, qseq_gens):
        """Yield qseq data"""
        _iter = qseq_gens
        for (seq_id, seq, qual, record) in _iter:
            self.state['SequenceID'] = seq_id
            self.state['Sequence'] = seq
            self.state['QualID'] = seq_id
            self.state['Qual'] = qual

            # as we're updating state in place and effectively circumventing
            # Workflow.initialize_state, we do not need to yield anything
            yield None
