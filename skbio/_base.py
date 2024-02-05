# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import abc


class SkbioObject(metaclass=abc.ABCMeta):
    """Abstract base class defining core API common to all scikit-bio objects.

    Public scikit-bio classes should subclass this class to ensure a common,
    core API is present. All abstract methods and properties defined here must
    be implemented in subclasses, otherwise they will not be instantiable.

    """

    @abc.abstractmethod
    def __str__(self):
        raise NotImplementedError


class ElasticLines:
    """Store blocks of content separated by dashed lines.

    Each dashed line (separator) is as long as the longest content
    (non-separator) line.

    """

    def __init__(self):
        self._lines = []
        self._separator_idxs = []
        self._max_line_len = -1

    def add_line(self, line):
        line_len = len(line)
        if line_len > self._max_line_len:
            self._max_line_len = line_len
        self._lines.append(line)

    def add_lines(self, lines):
        for line in lines:
            self.add_line(line)

    def add_separator(self):
        self._lines.append(None)
        self._separator_idxs.append(len(self._lines) - 1)

    def to_str(self):
        separator = "-" * self._max_line_len
        for idx in self._separator_idxs:
            self._lines[idx] = separator
        return "\n".join(self._lines)
