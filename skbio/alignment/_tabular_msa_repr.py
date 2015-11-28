# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from skbio.util._metadata_repr import _MetadataReprBuilder


class _TabularMSAReprBuilder(_MetadataReprBuilder):

    def _process_header(self):
        cls_name = self._obj.__class__.__name__
        dtype_class_name = self._obj._seqs[0].__class__.__name__
        self._lines.add_line(cls_name + ' <' + dtype_class_name + '>')
        self._lines.add_separator()

    def _process_data(self):
        num_lines = self._obj.shape.sequence

        # display entire sequence if we can, else display the first two and
        # last two lines separated by ellipsis
        if num_lines <= 5:
            self._lines.add_lines(self._format_lines(range(num_lines)))
        else:
            self._lines.add_lines(self._format_lines(range(2)))
            self._lines.add_line('...')
            self._lines.add_lines(self._format_lines(range(num_lines - 2,
                                                     num_lines)))

    def _format_lines(self, line_indices):
        num_chars_at_end = 4
        lines = []
        for line_index in line_indices:
            seq_str = str(self._obj._seqs[line_index])
            formatted_seq = (seq_str[0:-num_chars_at_end] + ' ... ' +
                             seq_str[-num_chars_at_end:])
            lines.append(formatted_seq)
        return lines
