#!/usr/bin/env python
r"""
Sequence collections and alignments (:mod:`bipy.core.alignment`)
================================================================

.. currentmodule:: bipy.core.alignment

This module provides functionality for working with biological sequence
collections and alignments. These can be composed of generic sequences, 
nucelotide sequences, DNA sequences, and RNA sequences.

Classes
-------

.. autosummary::
   :toctree: generated/

   SequenceCollection
   Alignment

Examples
--------
>>> from bipy.core.alignment import SequenceCollection, Alignment
>>> from bipy.core.sequence import DNA

>>> seqs = [("s1", "ACC--G-GGTA.."), ("s2", "TCC--G-GGCA..")]
>>> a1 = Alignment(seqs, DNA)

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


class SequenceCollection(object):
    """
    """

    def __init__(self, seqs, seq_constructor):
        """
        """
        self._data = []
        for seq_id, seq in seqs:
            try:
                identifier, description = seq_id.split(None, 1)
            except ValueError:
                identifier = seq_id.split()
                description = None
            self._data.append(seq_constructor(seq, identifier=identifier,
              description=description))

    def to_fasta(self):
        """
        """
        return ''.join([seq.to_fasta() for seq in self._data])


Alignment = SequenceCollection

