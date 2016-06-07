# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import math

from skbio.util._misc import chunk_str
from skbio.metadata._repr import _MetadataReprBuilder


class _SequenceReprBuilder(_MetadataReprBuilder):
    """Build a ``Sequence`` repr.

    Parameters
    ----------
    seq : Sequence
        Sequence to repr.
    width : int
        Maximum width of the repr.
    indent : int
        Number of spaces to use for indented lines.
    chunk_size: int
        Number of characters in each chunk of a sequence.

    """
    def __init__(self, seq, width, indent, chunk_size):
        super(_SequenceReprBuilder, self).__init__(seq, width, indent)
        self._chunk_size = chunk_size

    def _process_header(self):
        cls_name = self._obj.__class__.__name__
        self._lines.add_line(cls_name)
        self._lines.add_separator()

    def _process_data(self):
        num_lines, num_chars, column_width = self._find_optimal_seq_chunking()

        # display entire sequence if we can, else display the first two and
        # last two lines separated by ellipsis
        if num_lines <= 5:
            self._lines.add_lines(self._format_chunked_seq(
                range(num_lines), num_chars, column_width))
        else:
            self._lines.add_lines(self._format_chunked_seq(
                range(2), num_chars, column_width))
            self._lines.add_line('...')
            self._lines.add_lines(self._format_chunked_seq(
                range(num_lines - 2, num_lines), num_chars, column_width))

    def _find_optimal_seq_chunking(self):
        """Find the optimal number of sequence chunks to fit on a single line.

        Returns the number of lines the sequence will occupy, the number of
        sequence characters displayed on each line, and the column width
        necessary to display position info using the optimal number of sequence
        chunks.

        """
        # strategy: use an iterative approach to find the optimal number of
        # sequence chunks per line. start with a single chunk and increase
        # until the max line width is exceeded. when this happens, the previous
        # number of chunks is optimal
        num_lines = 0
        num_chars = 0
        column_width = 0

        num_chunks = 1
        not_exceeded = True
        while not_exceeded:
            line_len, new_chunk_info = self._compute_chunked_seq_line_len(
                num_chunks)
            not_exceeded = line_len <= self._width
            if not_exceeded:
                num_lines, num_chars, column_width = new_chunk_info
                num_chunks += 1
        return num_lines, num_chars, column_width

    def _compute_chunked_seq_line_len(self, num_chunks):
        """Compute line length based on a number of chunks."""
        num_chars = num_chunks * self._chunk_size

        # ceil to account for partial line
        num_lines = int(math.ceil(len(self._obj) / num_chars))

        # position column width is fixed width, based on the number of
        # characters necessary to display the position of the final line (all
        # previous positions will be left justified using this width)
        column_width = len('%d ' % ((num_lines - 1) * num_chars))

        # column width + number of sequence characters + spaces between chunks
        line_len = column_width + num_chars + (num_chunks - 1)
        return line_len, (num_lines, num_chars, column_width)

    def _format_chunked_seq(self, line_idxs, num_chars, column_width):
        """Format specified lines of chunked sequence data."""
        lines = []
        for line_idx in line_idxs:
            seq_idx = line_idx * num_chars
            chars = str(self._obj[seq_idx:seq_idx+num_chars])
            chunked_chars = chunk_str(chars, self._chunk_size, ' ')
            lines.append(('%d' % seq_idx).ljust(column_width) + chunked_chars)
        return lines
