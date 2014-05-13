#!/usr/bin/env python
r"""
Sequence factory (:mod:`skbio.factory.sequence`)
================================================

.. currentmodule:: skbio.factory.sequence

The sequence factory provides a standard mechanism to iterate over sequence
files regardless of file type or whether they are compressed. The method
will take a single or list of file paths, resolve the necessary file openers,
necessary iterator objects and the corresponding parsers.

Methods
-------

.. autosummary::
   :toctree: generated/

   factory

See Also
--------

skbio.core.iterator.SequenceIterator

Examples
--------

Since the ``factory`` is intended to operate on file paths, lets create two
files for the ``factory`` to use. The first one will be a regular FASTA file,
and the second will be a gzip'd FASTQ file:

>>> import os
>>> import gzip
>>> out = open('test_seqs.fna', 'w')
>>> out.write(">s1\nATGC\n>s2\nATGGC\n")
>>> out.close()
>>> outgz = gzip.open('test_seqs.fq.gz', 'w')
>>> _ = outgz.write("@s3\nAATTGG\n+\nghghgh\n@s4\nAAA\n+\nfgh\n")
>>> outgz.close()

Now, lets see what the factory can do:

>>> from skbio.factory.sequence import factory
>>> it = factory(['test_seqs.fna', 'test_seqs.fq.gz'], phred_offset=64)
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

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
from itertools import chain
from gzip import open as gzip_open

from future.utils.six.moves import zip

from skbio.core.iterator import FastaIterator, FastqIterator


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


def factory(seqs, qual=None, constructor=None, **kwargs):
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
    SequenceIterator instance
        the return is ``Iterable``

    See Also
    --------
    skbio.core.iterator.SequenceIterator
    skbio.core.iterator.FastaIterator
    skbio.core.iterator.FastqIterator
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
        raise ValueError("Cannot resolve iterators!")

    if _is_single_iterator_type(i_seqs):
        seqs_constructor = i_seqs[0]
        gen = seqs_constructor(seq=seqs, qual=qual, **kwargs)
    else:
        gen = chain(*[c(seq=[fp], **kwargs) for c, fp in zip(i_seqs, seqs)])

    return gen
