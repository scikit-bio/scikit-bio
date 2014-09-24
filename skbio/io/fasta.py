"""
FASTA format (:mod:`skbio.io.fasta`)
====================================

.. currentmodule:: skbio.io.fasta

TODO add description

Format Specification
--------------------
TODO add format specification

Format Parameters
-----------------
TODO add format parameters

Examples
--------
TODO add examples

References
----------
TODO add references

http://en.wikipedia.org/wiki/FASTA_format
http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.io import (register_reader, register_writer, register_sniffer,
                      FASTAFormatError)
from skbio.io._base import _chunk_str


@register_sniffer('fasta')
def _fasta_sniffer(obj, fh):
    pass


@register_reader('fasta')
def _fasta_to_generator(obj, fh):
    pass


@register_writer('fasta')
def _generator_to_fasta(obj, fh, max_width=None):
    for idx, seq in enumerate(obj):
        if len(seq) < 1:
            raise FASTAFormatError(
                "Cannot write biological sequence number %d in FASTA format "
                "because it does not contain any characters (i.e., it is an "
                "empty/blank sequence). Empty sequences are not supported in "
                "the FASTA file format." % (idx + 1))

        if seq.description:
            header = '%s %s' % (seq.id, seq.description)
        else:
            header = seq.id

        seq_str = str(seq)

        if max_width is not None:
            seq_str = _chunk_str(seq_str, max_width, '\n')

        fh.write('>%s\n%s\n' % (header, seq_str))
