#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os.path
from itertools import chain, izip
from numpy import array

from bipy.parse.fasta import MinimalFastaParser

### these need to move to bipy
from cogent.parse.fastq import MinimalFastqParser
from qiime.parse import is_casava_v180_or_later, MinimalQualParser
from qiime.quality import ascii_to_phred33, ascii_to_phred64
from qiime.util import gzip_open
###

from .workflow import Workflow, not_none

#### TODO:
#### if chain is reimplemented, it would be very easy to track state
#### on which file-like object contains a bad seq, what record, etc.

from collections import namedtuple

class SequenceRecord(namedtuple("SequenceRecord", ("SequenceID", "Sequence",
                                                   "QualID", "Qual"))):
    __slots__ = ()

    def __new__(cls, SequenceID, Sequence, QualID=None, Qual=None):
        return super(SequenceRecord, cls).__new__(cls, SequenceID, Sequence,
                                                  QualID, Qual)


### should these be member methods?
def valid_ids(a, b):
    """Check if a and b are equal"""
    return a == b


def valid_length(a, b):
    """Check if a and b are equal length"""
    return len(a) == len(b)


def has_qual(item):
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
    def __init__(self, seq, qual=None, transform=None, check_fail=None,
                 valid_id=True, valid_length=True):
        self.seq = seq
        self.qual = qual

        self._transform = transform
        self._check_fail = check_fail

        state = {'SequenceID': None,
                 'Sequence': None,
                 'QualID': None,
                 'Qual': None}

        options = {'transform': self._transform,
                   'check_fail': self._check_fail,
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

    @Workflow.method(priority=100)
    @Workflow.requires(option='valid_id', values=True)
    def validate_ids(self):
        self.valid_qual_id()

    @Workflow.method(priority=90)
    @Workflow.requires(option='valid_length', values=True, state=has_qual)
    def valid_lengths(self):
        self.failed = not valid_length(self.state['Sequence'],
                                       self.state['Qual'])

    @Workflow.method(priority=80)
    @Workflow.requires(option='transform', values=not_none)
    def transform(self):
        """Transform state if necessary"""
        self._transform(self.state)

    ## not sure if necessary
    @Workflow.method(priority=70)
    @Workflow.requires(option='check_fail', values=True)
    def check_fail(self):
        self.failed = self._check_fail(self.state)

    @Workflow.requires(state=has_qual)
    def valid_qual_id(self):
        self.failed = not valid_ids(self.state['SequenceID'],
                                    self.state['QualID'])

class FastaIterator(SequenceIterator):
    """Populate and yield records based on fasta sequence"""
    def _gen(self):
        """Construct internal iterators"""
        # construct fasta generators
        fasta_gens = chain(*[MinimalFastaParser(f) for f in self.seq])

        # construct qual generators if necessary
        if self.qual is not None:
            qual_gens = chain(*[MinimalQualParser(f) for f in self.qual])
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
        fastq_gens = chain(*[MinimalFastqParser(f, False) for f in self.seq])

        # peek the first record to determine phred offset
        first_item = fastq_gens.next()
        seqid, seq, qual = first_item
        fastq_gens = chain([first_item], fastq_gens)

        # from qiime.parse.parse_fastq_qual_score (v1.8.0)
        if is_casava_v180_or_later('@%s' % seqid):
            ascii_to_phred_f = ascii_to_phred33
        else:
            ascii_to_phred_f = ascii_to_phred64

        gen = self._fastq_gen(fastq_gens, ascii_to_phred_f)

        return gen

    def _fastq_gen(self, fastq_gens, phred_f):
        """Yield fastq data"""
        for (seq_id, seq, qual_str) in fastq_gens:
            qual = array([phred_f(q) for q in qual_str])
            yield SequenceRecord(seq_id, seq, seq_id, qual)

##### UNTESTED STUBBED OUT CODE BELOW #####

class PairedSequenceIterator(object):
    """Support for paired data from two sources

    The IDs of read N from the first source and read N from the second
    source are checked to make sure they are identical
    """
    def __init__(self, read1_iter, read2_iter):
        self.read1_iter = read1_iter
        self.read2_iter = read2_iter

    def __call__(self):
        for read1, read2 in izip(self.read1_iter, self.read2_iter):
            # is any unmangling of IDs necessary here?
            if not valid_ids(read1['SequenceID'], read2['SequenceID']):
                raise ValueError("%s and %s don't look paired!" % \
                        (read1['SequenceID'], read2['SequenceID']))

            yield (read1, read2)

class PairedShuffledSequenceIterator(object):
    """Support for paired read data that is shuffled

    The IDs of N and N+1 are checked to make sure they are identical
    """
    def __init__(self, reads_iter):
        self.reads_iter = reads_iter
        self._last_read = reads_iter.next()
        self.counter = 1

    # this likely needs to be abstracted further to handle other types of
    # abuses
    def _chop_pair_number(self, id_):
        """Removes the trailing "/0" or "/1" from a sequence ID"""
        return id_.rsplit('/', 1)[0]

    def __call__(self):
        last_read = self._last_read
        counter = self.counter

        for read in self.reads_iter:
            counter += 1

            if counter % 2 == 0:
                id1 = self._chop_pair_number(last_read['SequenceID'])
                id2 = self._chop_pair_number(read['SequenceID'])
                valid_ids(id1, id2)
                yield (last_read, read)
            else:
                last_read = read

        if counter % 2 != 0:
            raise ValueError("Uneven number of sequences observed!")

FILEEXT_MAP = {'fna':(FastaIterator, open),
               'fna.gz':(FastaIterator, gzip_open),
               'fasta':(FastaIterator, open),
               'fasta.gz':(FastaIterator, gzip_open),
               'qual':(FastaIterator, open),
               'qual.gz':(FastaIterator, gzip_open),
               'fastq':(FastqIterator, open),
               'fastq.gz':(FastqIterator, gzip_open)}

def _determine_types_and_openers(files):
    """Attempt to determine the appropriate iterators and openers"""
    if files is None:
        return [], []

    iters = []
    openers = []
    for fpath in files:
        if fpath.endswith('.gz'):
            ext = fpath.rsplit('.', 2)[-1]
        else:
            ext = fpath.rsplit('.', 1)[-1]

        i, o = FILEEXT_MAP.get(ext, None)
        if i is None:
            raise IOError("Unknown filetype for %s" % fpath)

        iters.append(i)
        openers.append(o)

    return iters, openers

def _is_single_iterator_type(iters):
    """Determine if there is a single or multiple type of iterator

    If iters is [], this method returns True it considers the null case to be
    a single iterator type.
    """
    if iters:
        return len(set(iters)) == 1
    else:
        return True

def _open_or_none(openers, files):
    """Open all the files or return None"""
    if not openers:
        return None
    else:
        opened = []
        for o, f in zip(openers, files):
            name = o.__name__
            if not os.path.exists(f):
                raise IOError("%s does not appear to exist!" % f)

            try:
                opened.append(o(f))
            except IOError as e:
                raise IOError("Could not open %s with %s!" % (f, name))

        return opened

def iter_factory(seqs1=None, qual1=None, bc1=None,
                 seqs2=None, qual2=None, bc2=None,
                 paired_and_shuffled=False, r1_kw=None, r2_kw=None):
    """Construct the appropriate iterator for all your processing needs

    This method will attempt to open all files correctly and to feed the
    appropriate objects into the correct iterators.

    seqs1 : list of sequence file paths for unpaired data or the first read if
        the data are paired
    seqs2 : list of sequence file paths for the second read of paired data
    qual1 : list of qual file paths for unpaired data or the first read if the
        data are paired
    qual2 : list of qual file paths for the second read of paired data
    bc1   : list of barcode file paths for unpaired data or the first read if
        the data are paired
    bc2   : list of barcode file paths for the second read of paired data
    paired_and_shuffled : if True, the factory will assume the data in seqs1,
        qual1 (if present) and bc1 (if present) are paired and interleaved
        internally
    r1_kw : a dict of arguments for the read 1 iterator, including transforms
        (e.g., reverse complement), validation (e.g., disable ID and length
        checks), etc.
    r2_kw : a dict of arguments for the read 2 iterator, including transforms
        (e.g., reverse complement), validation (e.g., disable ID and length
        checks), etc.

    Currently, the iterators cannot handle having a file as fastq and a file as
    fasta. In otherwords, the following will raise ValueError:

        x = iter_factory(seqs1=['foo.fasta', 'bar.fastq'])

    It is valid on paired data though, where the first read can be entirely
    fasta and the second can be entirely fastq:

        x = iter_factory(seqs1=['foo1.fasta', 'bar1.fasta'],
                         seqs2=['foo2.fastq', 'bar2.fastq'])
    """
    if seqs1 is None:
        raise ValueError("Must pass in sequences!")

    if r1_kw is None:
        r1_kw = {}

    if r2_kw is None:
        r2_kw = {}

    # i -> iters, o -> openers
    i_seqs1, o_seqs1 = _determine_types_and_openers(seqs1)
    i_seqs2, o_seqs2 = _determine_types_and_openers(seqs2)
    i_qual1, o_qual1 = _determine_types_and_openers(qual1)
    i_qual2, o_qual2 = _determine_types_and_openers(qual2)
    i_bc1, o_bc1 = _determine_types_and_openers(bc1)
    i_bc2, o_bc2 = _determine_types_and_openers(bc2)

    if seqs1 and not _is_single_iterator_type(i_seqs1):
        raise ValueError("Cannot mix and match iterators!")

    if seqs2 and not _is_single_iterator_type(i_seqs2):
        raise ValueError("Cannot mix and match iterators!")

    if qual1 and not _is_single_iterator_type(i_qual1):
        raise ValueError("Cannot mix and match iterators!")

    if qual2 and not _is_single_iterator_type(i_qual2):
        raise ValueError("Cannot mix and match iterators!")

    if bc1 and not _is_single_iterator_type(i_bc1):
        raise ValueError("Cannot mix and match iterators!")

    if bc2 and not _is_single_iterator_type(i_bc2):
        raise ValueError("Cannot mix and match iterators!")

    # fetch the constructor, and open all the files
    r1_seqs_constructor = i_seqs1[0]
    r1_seqs = _open_or_none(seqs1, o_seqs1)
    r1_qual = _open_or_none(qual1, o_qual1)
    r1_bc = _open_or_none(bc1, o_bc1)

    r2_seqs_constructor = i_seqs2[0] if i_seqs2 else None
    r2_seqs = _open_or_none(seqs2, o_seqs2)
    r2_qual = _open_or_none(qual2, o_qual2)
    r2_bc = _open_or_none(bc2, o_bc2)

    # build the read1 iterator
    r1 = r1_seqs_constructor(seqs=r1_seqs, qual=r1_qual, bc=r1_bc, **r1_kw)

    # handle paired reads if necessary
    if r2_seqs_constructor is not None:
        r2 = r2_seqs_constructor(seqs=r2_seqs, qual=r2_qual, bc=r2_bc, **r2_kw)
        obj = PairedSequenceIterator(r1, r2)
    elif paired_and_shuffled:
        obj = PairedShuffledSequenceIterator(r1)
    else:
        obj = r1

    return obj
