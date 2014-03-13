#!/usr/bin/env python
"""
Minimal Fasta Parser (:mod:`bipy.parse.fastq`)
=============================================

.. currentmodule:: bipy.fastq.parse

This module provides a function for parsing a fastq file.

Functions
---------

.. autosummary::
   :toctree: generated/

   MinimalFastqParser

"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


def MinimalFastqParser(data, strict=True):
    """yields name, seq, qual from fastq file

    Arguments:
        - strict: checks the quality and sequence labels are the same
    """

    # fastq format is very simple, defined by blocks of 4 lines
    line_num = -1
    record = []
    for line in data:
        line_num += 1
        if line_num == 4:
            if strict:  # make sure the seq and qual labels match
                assert record[0][1:] == record[2][1:], \
                    'Invalid format: %s -- %s' % (record[0][1:], record[2][1:])
            yield record[0][1:], record[1], record[3]

            line_num = 0
            record = []

        record.append(line.strip())

    if record:
        if strict and record[0]:  # make sure the seq and qual labels match
            assert record[0][1:] == record[2][1:], 'Invalid format'

        if record[0]:  # could be just an empty line at eof
            yield record[0][1:], record[1], record[3]

    if type(data) == file:
        data.close()
