from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from skbio.io import FileFormatError


class FastqParseError(FileFormatError):
    """Exception raised when a FASTQ formatted file cannot be parsed"""
    pass


class QseqParseError(FileFormatError):
    """Exception raised when a QSEQ formatted file cannot be parsed"""
    pass
