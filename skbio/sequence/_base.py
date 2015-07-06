# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from skbio.util._decorator import stable


class ElasticLines(object):
    """Store blocks of content separated by dashed lines.

    Each dashed line (separator) is as long as the longest content
    (non-separator) line.

    """

    @stable(as_of="0.4.0")
    def __init__(self):
        self._lines = []
        self._separator_idxs = []
        self._max_line_len = -1

    @stable(as_of="0.4.0")
    def add_line(self, line):
        line_len = len(line)
        if line_len > self._max_line_len:
            self._max_line_len = line_len
        self._lines.append(line)

    def add_lines(self, lines):
        for line in lines:
            self.add_line(line)

    @stable(as_of="0.4.0")
    def add_separator(self):
        self._lines.append(None)
        self._separator_idxs.append(len(self._lines) - 1)

    @stable(as_of="0.4.0")
    def to_str(self):
        separator = '-' * self._max_line_len
        for idx in self._separator_idxs:
            self._lines[idx] = separator
        return '\n'.join(self._lines)
