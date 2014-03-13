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

from bipy.core.exception import SequenceCollectionError

class SequenceCollection(object):
    """
    """

    def __init__(self, seqs, seq_constructor, validate=False):
        """
        """
        self._data = []
        for seq_id, seq in seqs:
            try:
                identifier, description = seq_id.split(None, 1)
            except ValueError:
                identifier = seq_id.strip()
                description = None
            self._data.append(seq_constructor(seq, identifier=identifier,
              description=description))
 
        if validate and not self.is_valid():
            raise SequenceCollectionError(
                "Something is wrong, and it's your fault.")
    
    def __getitem__(self, index):
        """
        """
        return self._data[index]

    def __iter__(self):
        """
        """
        return iter(self._data)

    def __len__(self):
        """
        """
        return len(self._data)

    def _validate_character_set(self):
        """
        """
        for seq in self:
            if not seq.is_valid():
                return False
        return True

    def is_valid(self):
        """
        """
        return self._validate_character_set()

    def to_fasta(self):
        """
        """
        return ''.join([seq.to_fasta() for seq in self._data])


class Alignment(SequenceCollection):
    """
    """

    def _validate_lengths(self):
        """
        """
        seq1_length = len(self[0])
        for seq in self:
            if seq1_length != len(seq):
                return False
        return True

    def is_valid(self):
        """
        """
        result = super(Alignment, self).is_valid() or self._validate_lengths()

