#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import os
from gzip import open as gzip_open
from .iterator import FastaIterator, FastqIterator

FILEEXT_MAP = {'fna': (FastaIterator, open),
               'fna.gz': (FastaIterator, gzip_open),
               'fasta': (FastaIterator, open),
               'fasta.gz': (FastaIterator, gzip_open),
               'qual': (FastaIterator, open),
               'qual.gz': (FastaIterator, gzip_open),
               'fastq': (FastqIterator, open),
               'fastq.gz': (FastqIterator, gzip_open)}

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


def factory(seqs, qual=None, **kwargs):
    """Construct the appropriate iterator for all your processing needs

    This method will attempt to open all files correctly and to feed the
    appropriate objects into the correct iterators.

    seqs : list of sequence file paths
    qual : list of qual file paths or None

    kwargs are passed into the subsequent generators.

    Seqs can list multiple types of files (e.g., FASTA and FASTQ), but if
    multiple file types are specified, qual must be None
    """
    if not seqs:
        raise ValueError("Must pass in sequences!")

    # i -> iters, o -> openers
    i_seqs, o_seqs = _determine_types_and_openers(seqs)
    i_qual, o_qual = _determine_types_and_openers(qual)

    seqs = [_open_or_none(f) for f in seqs]
    qual = [_open_or_none(f) for f in qual] or None

    if not _is_single_iterator_type(i_seqs) and qual is not None:
        # chaining Fasta/Fastq for sequence is easy, but it gets nasty quick
        # if seqs is a mix of fasta/fastq, with qual coming in as there aren't
        # 1-1 mappings. This could be addressed if necessary, but seems like
        # an unnecessary block of code right now
        raise ValueError("Cannot resolve iterators!")

    if _is_single_iterator_type(i_seqs):
        seqs_constructor = i_seqs[0]
        gen = seqs_constructor(seqs=seqs, qual=qual, **kwargs)
    else:
        gen = chain(*[c(seqs=fp, **kwargs) for c, fp in izip(i_seqs, seqs)])

    return gen
