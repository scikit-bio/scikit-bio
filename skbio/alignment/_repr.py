# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util._metadata_repr import _MetadataReprBuilder
import math


class _TabularMSAReprBuilder(_MetadataReprBuilder):

    def _process_header(self):
        cls_name = self._obj.__class__.__name__
        if self._obj.dtype is not None:
            dtype_class = '<' + self._obj.dtype.__name__ + '>'
        else:
            dtype_class = ''
        self._lines.add_line(cls_name + dtype_class)
        self._lines.add_separator()
        self._ellipse_insert = ' ... '

    def _process_data(self):
        num_sequences = self._obj.shape.sequence
        num_positions = self._obj.shape.position

        # catch case of all empty sequences
        if num_positions > 0:
            # display all sequences if we can, else display the first two and
            # last two sequences separated by ellipsis
            if num_sequences <= 5:
                self._lines.add_lines(
                    self._format_sequences(range(num_sequences)))
            else:
                self._lines.add_lines(self._format_sequences(range(2)))
                self._lines.add_line('...')
                self._lines.add_lines(self._format_sequences(
                    range(num_sequences - 2, num_sequences)))

    def _format_sequences(self, line_indices):
        lines = []
        for line_index in line_indices:
            seq_str = str(self._obj._get_sequence(line_index))
            if len(seq_str) <= self._width:
                formatted_seq = seq_str
            else:
                formatted_seq = (
                    seq_str[0:self._num_characters_before_ellipse()] +
                    self._ellipse_insert +
                    seq_str[-self._num_characters_after_ellipse():]
                )
            lines.append(formatted_seq)
        return lines

    def _num_characters_before_ellipse(self):
        return math.floor(self._num_characters_to_display() / 2)

    def _num_characters_after_ellipse(self):
        return (self._num_characters_to_display() -
                self._num_characters_before_ellipse())

    def _num_characters_to_display(self):
        return self._width - len(self._ellipse_insert)