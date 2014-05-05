r"""
Format biological sequences (:mod:`skbio.format.sequences`)
===========================================================

.. currentmodule:: skbio.format.sequences

This module provides functions for writing sequence files in a variety of
different formats, the available formatters are listed below.

Functions
---------

.. autosummary::
   :toctree: generated/

    fasta_from_sequences
    fasta_from_alignment
    format_fastq_record'

"""
#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from .fasta import fasta_from_sequences, fasta_from_alignment
from .fastq import format_fastq_record

__all__ = ['fasta_from_sequences', 'fasta_from_alignment',
           'format_fastq_record']


from numpy.testing import Tester
test = Tester().test
