#!/usr/bin/env python
r"""
Parse biological sequences (:mod:`skbio.parse.sequences`)
=========================================================

.. currentmodule:: skbio.parse.sequences

This module provides functions for parsing sequence files in a variety of
different formats, the available parsers are listed below.

Functions
---------

.. autosummary::
   :toctree: generated/

    parse_clustal
    parse_fasta
    parse_fastq


"""
# -----------------------------------------------------------------------------
# Copyright (c) 2014--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from .clustal import parse_clustal
from .fasta import parse_fasta
from .fastq import parse_fastq

__all__ = ['parse_clustal', 'parse_fasta', 'parse_fastq']
