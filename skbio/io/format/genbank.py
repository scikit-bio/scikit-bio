"""
Genbank format (:mod:`skbio.io.format.genbank`)
===============================================

.. currentmodule:: skbio.io.format.genbank

Format Support
--------------
**Has Sniffer: Yes**

+------+------+---------------------------------------------------------------+
|Reader|Writer|                          Object Class                         |
+======+======+===============================================================+
|Yes   |Yes   |:mod:`skbio.sequence.Sequence`                                 |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.DNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.RNA`                                      |
+------+------+---------------------------------------------------------------+
|Yes   |Yes   |:mod:`skbio.sequence.Protein`                                  |
+------+------+---------------------------------------------------------------+

Format Specification
--------------------
The following sections define the genbank formats.

<http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>

Feature Table Documentation:
http://www.insdc.org/files/feature_table.html
http://www.ncbi.nlm.nih.gov/projects/collab/FT/index.html
ftp://ftp.ncbi.nih.gov/genbank/docs/

    The International Nucleotide Sequence Database Collaboration (INSDC)
    between the DDBJ, EMBL, and GenBank.  These organisations all use the
    same "Feature Table" layout in their plain text flat file formats.

    However, the header and sequence sections of an EMBL file are very
    different in layout to those produced by GenBank/DDBJ.

Examples
--------

"""
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import range, zip
from six.moves import zip_longest
import re
import textwrap

import numpy as np

from skbio.io import create_format, GenbankFormatError
from skbio.io.registry import FileSentinel
from skbio.io.format._base import (_get_nth_sequence,
                                   _parse_fasta_like_header,
                                   _format_fasta_like_records, _line_generator,
                                   _too_many_blanks)
from skbio.util._misc import chunk_str
from skbio.sequence import Sequence, DNA, RNA, Protein


genbank = create_format('genbank')


@genbank.sniffer()
def _genbank_sniffer(fh):
    # check the 1st line starts with 'LOCUS' and the last line ends with '//'
    if _too_many_blanks(fh, 5):
        return False, {}

    empty = True
    try:
        parser = _parse_genbank_raw(fh)
        for _ in parser:
            empty = False
    except GenbankFormatError:
        return False, {}
    if empty:
        return False, {}
    else:
        return True, {}


@genbank.reader(None)
def _genbank_to_generator(fh, constructor=Sequence, **kwargs):
    pass


@genbank.reader(Sequence)
def _genbank_to_biological_sequence(fh, rec_num=1):
    pass


@genbank.reader(DNA)
def _genbank_to_dna(fh, rec_num=1):
    pass


@genbank.reader(RNA)
def _genbank_to_rna(fh, rec_num=1):
    pass


@genbank.reader(Protein)
def _genbank_to_protein(fh, rec_num=1):
    pass


@genbank.writer(None)
def _sequences_to_genbank(fh):
    pass


def _parse_genbank_raw(fh):
    data_chunks = []
    for line in _line_generator(fh, skip_blanks=True):
        if line.startswith('//'):
            yield data_chunks
            data_chunks = []
        else:
            data_chunks.append(line)
