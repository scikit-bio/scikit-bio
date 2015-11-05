# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
import six

import math
import itertools
import numbers
import textwrap

from skbio.sequence._base import ElasticLines
from skbio.util._misc import chunk_str


class _SequenceReprBuilder(object):
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
        self._seq = seq
        self._width = width
        self._indent = ' ' * indent
        self._chunk_size = chunk_size

    def build(self):
        lines = ElasticLines()

        cls_name = self._seq.__class__.__name__
        lines.add_line(cls_name)
        lines.add_separator()

        if self._seq.has_metadata():
            lines.add_line('Metadata:')
            # Python 3 doesn't allow sorting of mixed types so we can't just
            # use sorted() on the metadata keys. Sort first by type then sort
            # by value within each type.
            for key in self._sorted_keys_grouped_by_type(self._seq.metadata):
                value = self._seq.metadata[key]
                lines.add_lines(self._format_metadata_key_value(key, value))

        if self._seq.has_positional_metadata():
            lines.add_line('Positional metadata:')
            for key in self._seq.positional_metadata.columns.values.tolist():
                dtype = self._seq.positional_metadata[key].dtype
                lines.add_lines(
                    self._format_positional_metadata_column(key, dtype))

        lines.add_line('Stats:')
        for label, value in self._seq._repr_stats():
            lines.add_line('%s%s: %s' % (self._indent, label, value))
        lines.add_separator()

        num_lines, num_chars, column_width = self._find_optimal_seq_chunking()

        # display entire sequence if we can, else display the first two and
        # last two lines separated by ellipsis
        if num_lines <= 5:
            lines.add_lines(self._format_chunked_seq(
                range(num_lines), num_chars, column_width))
        else:
            lines.add_lines(self._format_chunked_seq(
                range(2), num_chars, column_width))
            lines.add_line('...')
            lines.add_lines(self._format_chunked_seq(
                range(num_lines - 2, num_lines), num_chars, column_width))

        return lines.to_str()

    def _sorted_keys_grouped_by_type(self, dict_):
        """Group keys within a dict by their type and sort within type."""
        type_sorted = sorted(dict_, key=self._type_sort_key)
        type_and_value_sorted = []
        for _, group in itertools.groupby(type_sorted, self._type_sort_key):
            type_and_value_sorted.extend(sorted(group))
        return type_and_value_sorted

    def _type_sort_key(self, key):
        return repr(type(key))

    def _format_metadata_key_value(self, key, value):
        """Format metadata key:value, wrapping across lines if necessary."""
        key_fmt = self._format_key(key)

        supported_type = True
        if isinstance(value, (six.text_type, six.binary_type)):
            # for stringy values, there may be u'' or b'' depending on the type
            # of `value` and version of Python. find the starting quote
            # character so that wrapped text will line up with that instead of
            # the string literal prefix character. for example:
            #
            #     'foo': u'abc def ghi
            #              jkl mno'
            value_repr = repr(value)
            extra_indent = 1
            if not (value_repr.startswith("'") or value_repr.startswith('"')):
                extra_indent = 2
        # handles any number, this includes bool
        elif value is None or isinstance(value, numbers.Number):
            value_repr = repr(value)
            extra_indent = 0
        else:
            supported_type = False

        if not supported_type or len(value_repr) > 140:
            value_repr = str(type(value))
            # extra indent of 1 so that wrapped text lines up past the bracket:
            #
            #     'foo': <type
            #             'dict'>
            extra_indent = 1

        return self._wrap_text_with_indent(value_repr, key_fmt, extra_indent)

    def _format_key(self, key):
        """Format metadata key.

        Includes initial indent and trailing colon and space:

            <indent>'foo':<space>

        """
        key_fmt = self._indent + repr(key)
        supported_types = (six.text_type, six.binary_type, numbers.Number,
                           type(None))
        if len(key_fmt) > (self._width / 2) or not isinstance(key,
                                                              supported_types):
            key_fmt = self._indent + str(type(key))
        return '%s: ' % key_fmt

    def _wrap_text_with_indent(self, text, initial_text, extra_indent):
        """Wrap text across lines with an initial indentation.

        For example:

            'foo': 'abc def
                    ghi jkl
                    mno pqr'

        <indent>'foo':<space> is `initial_text`. `extra_indent` is 1. Wrapped
        lines are indented such that they line up with the start of the
        previous line of wrapped text.

        """
        return textwrap.wrap(
            text, width=self._width, expand_tabs=False,
            initial_indent=initial_text,
            subsequent_indent=' ' * (len(initial_text) + extra_indent))

    def _format_positional_metadata_column(self, key, dtype):
        key_fmt = self._format_key(key)
        dtype_fmt = '<dtype: %s>' % str(dtype)
        return self._wrap_text_with_indent(dtype_fmt, key_fmt, 1)

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
        num_lines = int(math.ceil(len(self._seq) / num_chars))

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
            chars = str(self._seq[seq_idx:seq_idx+num_chars])
            chunked_chars = chunk_str(chars, self._chunk_size, ' ')
            lines.append(('%d' % seq_idx).ljust(column_width) + chunked_chars)
        return lines
