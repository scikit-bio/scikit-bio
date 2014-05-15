r"""
Parse biological sequences (:mod:`skbio.parse.sequences`)
=========================================================

.. currentmodule:: skbio.parse.sequences

This module provides functions for parsing sequence files in a variety of
different formats, the available parsers are listed below.

This module also provides standardized iterators. The benefit of using these
iterators for things like sequence files is that the type of the file, and
any file format details are abstracted out to the developer. In this manner,
the developer does not need to worry about whether they're operating on
FASTA or FASTQ, and any differences in the returns from their respective
parsers.

The `load` function provides a standard mechanism to iterate over sequence
files regardless of file type or whether they are compressed. The method
will take a single or list of file paths, resolve the necessary file openers,
necessary iterator objects and the corresponding parsers.

Classes
-------

.. autosummary::
   :toctree: generated/

   SequenceIterator
   FastaIterator
   FastqIterator

Functions
---------

.. autosummary::
   :toctree: generated/

    parse_clustal
    parse_fasta
    parse_fastq
    parse_qual
    load

Examples
--------

For the first set of examples, we're going to use the `load` function. The
`load` function is intended to operate on file paths, so lets create two files
for it to use. The first one will be a regular FASTA file, and the second will
be a gzip'd FASTQ file:

>>> import os
>>> import gzip
>>> out = open('test_seqs.fna', 'w')
>>> out.write(">s1\nATGC\n>s2\nATGGC\n")
>>> out.close()
>>> outgz = gzip.open('test_seqs.fq.gz', 'w')
>>> _ = outgz.write("@s3\nAATTGG\n+\nghghgh\n@s4\nAAA\n+\nfgh\n")
>>> outgz.close()

Now, lets see what `load` can do:

>>> it = load(['test_seqs.fna', 'test_seqs.fq.gz'], phred_offset=64)
>>> for rec in it:
...     print rec['SequenceID']
...     print rec['Sequence']
...     print rec['Qual']
s1
ATGC
None
s2
ATGGC
None
s3
AATTGG
[39 40 39 40 39 40]
s4
AAA
[38 39 40]

To be polite, lets now remove the files we just created:

>>> os.remove('test_seqs.fna')
>>> os.remove('test_seqs.fq.gz')

The next set of examples surround the `SequenceIterator` and derived classes.

>>> from StringIO import StringIO
>>> from skbio.parse.sequences import FastaIterator, FastqIterator

In the first example, we're going to construct a FASTA iterator that is also
paired with quality scores (e.g., as in 454 fasta/qual files).

>>> seqs = StringIO(">seq1\n"
...                 "ATGC\n"
...                 ">seq2\n"
...                 "TTGGCC\n")
>>> qual = StringIO(">seq1\n"
...                 "10 20 30 40\n"
...                 ">seq2\n"
...                 "1 2 3 4 5 6\n")
>>> it = FastaIterator(seq=[seqs], qual=[qual])
>>> for record in it:
...     print record['Sequence']
...     print record['Qual']
ATGC
[10 20 30 40]
TTGGCC
[1 2 3 4 5 6]

In the next example, we're going to iterate over multiple FASTQ files at once.

>>> seqs1 = StringIO("@seq1\n"
...                  "ATGC\n"
...                  "+\n"
...                  "hhhh\n")
>>> seqs2 = StringIO("@seq2\n"
...                 "AATTGGCC\n"
...                 ">seq2\n"
...                 "abcdefgh\n")
>>> it = FastqIterator(seq=[seqs1, seqs2], phred_offset=64)
>>> for record in it:
...     print record['Sequence']
...     print record['Qual']
ATGC
[40 40 40 40]
AATTGGCC
[33 34 35 36 37 38 39 40]

Finally, we can apply arbitrary transforms to the sequences during iteratation.

>>> seqs1 = StringIO("@seq1\n"
...                  "ATGC\n"
...                  "+\n"
...                  "hhhh\n")
>>> seqs2 = StringIO("@seq2\n"
...                 "AATTGGCC\n"
...                 ">seq2\n"
...                 "abcdefgh\n")
>>> def rev_f(st):
...     st['Sequence'] = st['Sequence'][::-1]
...     st['Qual'] = st['Qual'][::-1] if st['Qual'] is not None else None
>>> it = FastqIterator(seq=[seqs1, seqs2], transform=rev_f, phred_offset=64)
>>> for record in it:
...     print record['Sequence']
...     print record['Qual']
CGTA
[40 40 40 40]
CCGGTTAA
[40 39 38 37 36 35 34 33]

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from itertools import chain
from gzip import open as gzip_open

from future.builtins import zip
from numpy.testing import Tester

from skbio.core.workflow import Workflow, not_none, method, requires

from .clustal import parse_clustal
from .fasta import parse_fasta, parse_qual
from .fastq import parse_fastq


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


def _determine_types_and_openers(files):
    """Attempt to determine the appropriate iterators and openers"""
    if files is None:
        return [], []

    iters = []
    openers = []
    for fpath in files:
        if fpath.endswith('.gz'):
            ext = '.'.join(fpath.rsplit('.', 2)[-2:])
        else:
            ext = fpath.rsplit('.', 1)[-1]

        i, o = FILEEXT_MAP.get(ext, (None, None))
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


def _open_or_none(opener, f):
    """Open a file or returns None"""
    if not opener:
        return None
    else:
        name = opener.__name__

        if not os.path.exists(f):
            raise IOError("%s does not appear to exist!" % f)

        try:
            opened = opener(f)
        except IOError:
            raise IOError("Could not open %s with %s!" % (f, name))

        return opened


def load(seqs, qual=None, constructor=None, **kwargs):
    """Construct the appropriate iterator for all your processing needs

    This method will attempt to open all files correctly and to feed the
    appropriate objects into the correct iterators.

    Seqs can list multiple types of files (e.g., FASTA and FASTQ), but if
    multiple file types are specified, qual must be None

    Parameters
    ----------
    seqs : str or list of sequence file paths
    qual : str or list of qual file paths or None
    constructor : force a constructor on seqs
    kwargs : dict
        passed into the subsequent generators.

    Returns
    -------
    SequenceIterator
        the return is ``Iterable``

    See Also
    --------
    skbio.parse.sequences.SequenceIterator
    skbio.parse.sequences.FastaIterator
    skbio.parse.sequences.FastqIterator
    """
    if not seqs:
        raise ValueError("Must pass in sequences!")

    if isinstance(seqs, str):
        seqs = [seqs]

    if isinstance(qual, str):
        qual = [qual]

    # i -> iters, o -> openers
    if constructor is not None:
        i_seqs = [constructor] * len(seqs)
        o_seqs = [open] * len(seqs)
    else:
        i_seqs, o_seqs = _determine_types_and_openers(seqs)

    i_qual, o_qual = _determine_types_and_openers(qual)

    seqs = [_open_or_none(o, f) for f, o in zip(seqs, o_seqs)]
    qual = [_open_or_none(o, f) for f, o in zip(qual or [], o_qual or [])]

    if not qual:
        qual = None

    if not _is_single_iterator_type(i_seqs) and qual is not None:
        # chaining Fasta/Fastq for sequence is easy, but it gets nasty quick
        # if seqs is a mix of fasta/fastq, with qual coming in as there aren't
        # 1-1 mappings. This could be addressed if necessary, but seems like
        # an unnecessary block of code right now
        raise ValueError("Cannot handle multiple sequence file types and qual "
                         "at the sametime!")

    if _is_single_iterator_type(i_seqs):
        seqs_constructor = i_seqs[0]
        gen = seqs_constructor(seq=seqs, qual=qual, **kwargs)
    else:
        gen = chain(*[c(seq=[fp], **kwargs) for c, fp in zip(i_seqs, seqs)])

    return gen


FILEEXT_MAP = {'fna': (FastaIterator, open),
               'fna.gz': (FastaIterator, gzip_open),
               'fasta': (FastaIterator, open),
               'fasta.gz': (FastaIterator, gzip_open),
               'qual': (FastaIterator, open),
               'qual.gz': (FastaIterator, gzip_open),
               'fastq': (FastqIterator, open),
               'fastq.gz': (FastqIterator, gzip_open),
               'fq': (FastqIterator, open),
               'fq.gz': (FastqIterator, gzip_open)}


__all__ = ['parse_clustal', 'parse_fasta', 'parse_fastq', 'parse_qual',
           'FastqIterator', 'FastaIterator', 'SequenceIterator', 'load']
